#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <byteswap.h>
#include <stdbool.h>

#include "jpeg-6b-steg/jpeglib.h"
#include "iterator.h"

struct jpeg_error_mgr jerr;

static void on_error_exit(j_common_ptr cinfo) {
	char buffer[JMSG_LENGTH_MAX];

	(*cinfo->err->format_message)(cinfo, buffer);

	jpeg_destroy(cinfo);

  puts(buffer);
}
static void on_output_message(j_common_ptr cinfo) {
	char buffer[JMSG_LENGTH_MAX];

	(*cinfo->err->format_message)(cinfo, buffer);

	puts(buffer);
}

static jvirt_barray_ptr* read_DCT_coefficients(struct jpeg_decompress_struct* cinfo, FILE* image) {
  cinfo->err = jpeg_std_error(&jerr);
  cinfo->err->error_exit = on_error_exit;
  cinfo->err->output_message = on_output_message;
  jpeg_create_decompress(cinfo);

  jpeg_stdio_src(cinfo, image);

  (void)jpeg_read_header(cinfo, TRUE);

  jpeg_calc_output_dimensions(cinfo);

  return jpeg_read_coefficients(cinfo);
}

typedef struct coeff_info {
  JCOEF* ptr; // Pointer to the coefficient in memory
  bool used;
  int target_lsb;
} coeff_info;

typedef struct bitmap {
  coeff_info* coeffs;
  int capacity;
  int count;
} bitmap;

bool bitmap_init(bitmap* bm) {
  memset(bm, 0, sizeof(*bm));

  bm->capacity = 1024;
  bm->coeffs = malloc(bm->capacity * sizeof(coeff_info));

  return bm->coeffs != NULL;
}
bool bitmap_append_bit(bitmap* bm, uint8_t bit, coeff_info origin) {
  bm->coeffs[bm->count++] = origin;

  if (bm->count < bm->capacity)
    return true;

  // Realloc
  bm->capacity += 1024;
  coeff_info* new_coeffs = realloc(bm->coeffs, bm->capacity * sizeof(coeff_info));

  if (new_coeffs == NULL) {
    free(new_coeffs);
    return false;
  }

  bm->coeffs = new_coeffs;
  return true;
}
void bitmap_free(bitmap* bm) {
  free(bm->coeffs);
}
void bitmap_reset_usage(bitmap* bm) {
  for (unsigned int i = 0; i < bm->count; i++)
    bm->coeffs[i].used = false;
}
coeff_info* bitmap_access_bit(bitmap* bm, unsigned int index) {
  return &bm->coeffs[index];
}

static int detect_quality_table(JQUANT_TBL* table, const unsigned int* basic_table) {
  for (int quality = 1; quality <= 100; quality++) {
    int percentage_scaling = quality < 50 ? 5000 / quality : 200 - quality*2;

    bool found_quality = true;
    for (int i = 0; i < DCTSIZE2; i++) {
      long expected_value = (basic_table[i] * percentage_scaling + 50) / 100;

      // Let's force baseline compatibility, as that's the default settings of JPEG 6b
      if (expected_value <= 0) expected_value = 1;
      if (expected_value > 255) expected_value = 255;
      
      if (table->quantval[i] != expected_value) {
        found_quality = false;
        break;
      }
    }

    if (found_quality)
      return quality;
  }

  return -1;
}
static int detect_quality(struct jpeg_decompress_struct* cinfo) {
  int quality = -1;

  // Luminance
  JQUANT_TBL* table = cinfo->quant_tbl_ptrs[0];
  if (table == NULL) {
    puts("Luminance table missing.");
    return -1;
  }

  static const unsigned int std_luminance_quant_tbl[DCTSIZE2] = {
    16,  11,  10,  16,  24,  40,  51,  61,
    12,  12,  14,  19,  26,  58,  60,  55,
    14,  13,  16,  24,  40,  57,  69,  56,
    14,  17,  22,  29,  51,  87,  80,  62,
    18,  22,  37,  56,  68, 109, 103,  77,
    24,  35,  55,  64,  81, 104, 113,  92,
    49,  64,  78,  87, 103, 121, 120, 101,
    72,  92,  95,  98, 112, 100, 103,  99
  };

  quality = detect_quality_table(table, std_luminance_quant_tbl);
  if (quality == -1) {
    puts("Could not detect quality of the luminance table.");
    return -1;
  }

  // Chrominance
  table = cinfo->quant_tbl_ptrs[1];
  if (table == NULL) {
    puts("Chrominance table missing.");
    return -1;
  }

  static const unsigned int std_chrominance_quant_tbl[DCTSIZE2] = {
    17,  18,  24,  47,  99,  99,  99,  99,
    18,  21,  26,  66,  99,  99,  99,  99,
    24,  26,  56,  99,  99,  99,  99,  99,
    47,  66,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99
  };

  int table_quality = detect_quality_table(table, std_chrominance_quant_tbl);
  if (table_quality != quality) {
    puts("Could not detect quality of the chrominance table.");
    return -1;
  }

  // Any other table?
  if (cinfo->quant_tbl_ptrs[2] != NULL) {
    puts("Extra quantization table, can't detect the quality of these.");
    return -1;
  }

  return quality;  
}

int min(int a, int b) {
  return a < b ? a : b;
}

// a length too small would result in a OOB read in Outguess... (:
#define DATA_LENGTH 128

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("Usage: %s <image.jpg> n\n", argv[0]);
    return 0;
  }

  char* input_filename = argv[1];
  FILE* input = fopen(input_filename, "rb");
  if (input == NULL) {
    printf("Cannot open '%s'.\n", input_filename);
    return 1;
  }

  /* Reading input image. */
  struct jpeg_decompress_struct in_cinfo;
  jvirt_barray_ptr* coeffs = read_DCT_coefficients(&in_cinfo, input);

  jpeg_component_info* components_info = in_cinfo.comp_info;

  // Storage for coefficients that can be used to hide data
  bitmap bm;
  if (!bitmap_init(&bm)) {
    puts("Out-of-memory error.");

    jpeg_finish_decompress(&in_cinfo);
    jpeg_destroy_decompress(&in_cinfo);
    fclose(input);

    return 1;
  }

  // Size in blocks of the downsampled Cb and Cr planes
  int width_blocks = components_info[1].width_in_blocks;
  int height_blocks = components_info[1].height_in_blocks;

  for (int by = 0; by < height_blocks; by++) {
    for (int bx = 0; bx < width_blocks; bx++) {
      jvirt_barray_ptr component_coeffs;
      JBLOCKROW row_coeffs;
      JCOEFPTR block_coeffs;

      for (int plane = 0; plane < 3; plane++) {
        int width_blocks = components_info[plane].width_in_blocks;
        int height_blocks = components_info[plane].height_in_blocks;
  
        int h_sampling_vector = components_info[plane].h_samp_factor;
        int v_sampling_vector = components_info[plane].v_samp_factor;
  
        for (int sub_by = v_sampling_vector*by; sub_by < min(height_blocks, v_sampling_vector*(by + 1)); sub_by++) {
          for (int sub_bx = h_sampling_vector*bx; sub_bx < min(width_blocks, h_sampling_vector*(bx + 1)); sub_bx++) {
            component_coeffs = coeffs[plane];
            row_coeffs = in_cinfo.mem->access_virt_barray((j_common_ptr)&in_cinfo, component_coeffs, sub_by, 1, TRUE)[0];
            block_coeffs = row_coeffs[sub_bx];
  
            for (int i = 0; i < DCTSIZE2; i++) {
              if (block_coeffs[i] == 0 || block_coeffs[i] == 1)
                continue;

              // Copy coefficient
              if (!bitmap_append_bit(&bm, block_coeffs[i], (coeff_info){ .used = false, .ptr = &block_coeffs[i] })) {
                puts("Out-of-memory error.");
  
                bitmap_free(&bm);
                jpeg_finish_decompress(&in_cinfo);
                jpeg_destroy_decompress(&in_cinfo);
                fclose(input);
  
                return 1;
              }
            }
          }
        }
      }
    }
  }

  printf("Usable coefficients: %d\n", bm.count);

  /* Find key prefix */
  int n = atoi(argv[2]);
  char* key_prefix = malloc(32);

  for (unsigned long long seed = 0; ; seed++) {
    sprintf(key_prefix, "%lld", seed);
    int key_prefix_length = strlen(key_prefix);
    key_prefix[key_prefix_length] = '_';
    key_prefix[key_prefix_length + 1] = '\0';
    key_prefix_length++;

    char* key = malloc(key_prefix_length + 32 + 1);
    strcpy(key, key_prefix);

    bitmap_reset_usage(&bm);
    bool is_key_suitable = true;

    for (int i = 0; i < n; i++) {
      sprintf(key + key_prefix_length, "%d", i);

      iterator it; // Iterator over the bitmap
      iterator_init(&it, key, strlen(key));

      struct arc4_stream as; // Encryption stream
      arc4_initkey(&as, "Encryption", key, strlen(key));

      coeff_info* header_bits[32];
      for (int j = 0; j < 32; j++) {
        header_bits[j] = bitmap_access_bit(&bm, ITERATOR_CURRENT(&it));
        iterator_next(&it);
      }

      bool is_ok = true;
      for (unsigned int j = 0; j < 4; j++) {
        unsigned char encryption_byte = arc4_getbyte(&as);

        // Don't care about the first two bytes: the seed.
        if (j < 2) continue;

        unsigned char data_byte = j == 2 ? DATA_LENGTH : 0;
        data_byte ^= encryption_byte;

        for (unsigned int k = 0; k < 8; k++) {
          unsigned int l = 8*j + k;

          int lsb = (data_byte >> k) & 1;

          if (header_bits[l]->used && header_bits[l]->target_lsb != lsb) {
            is_ok = false;
            break;
          }

          header_bits[l]->target_lsb = lsb;
          header_bits[l]->used = true;
        }
      }

      if (!is_ok) {
        is_key_suitable = false;
        break;
      }
    }

    if (is_key_suitable)
      break;

    free(key);
  }

  printf("Suitable key prefix: %s\n", key_prefix);

  /* Found suitable seed */
  int key_prefix_length = strlen(key_prefix);
  char* key = malloc(key_prefix_length + 32 + 1);
  strcpy(key, key_prefix);

  bitmap_reset_usage(&bm);
  bool is_key_suitable = true;

  for (int i = 0; i < n; i++) {
    sprintf(key + key_prefix_length, "%d", i);

    iterator it; // Iterator over the bitmap
    iterator_init(&it, key, strlen(key));

    struct arc4_stream as; // Encryption stream
    arc4_initkey(&as, "Encryption", key, strlen(key));

    coeff_info* header_bits[32];
    for (int j = 0; j < 32; j++) {
      header_bits[j] = bitmap_access_bit(&bm, ITERATOR_CURRENT(&it));
      iterator_next(&it);
    }

    bool is_ok = true;
    // Don't care about the first two bytes : the seed.
    for (unsigned int j = 0; j < 4; j++) {
      unsigned char encryption_byte = arc4_getbyte(&as);

      // Don't care about the first two bytes: the seed.
      if (j < 2) continue;

      unsigned char data_byte = j == 2 ? DATA_LENGTH : 0;
      data_byte ^= encryption_byte;

      for (unsigned int k = 0; k < 8; k++) {
        unsigned int l = 8*j + k;

        int lsb = (data_byte >> k) & 1;
        *(header_bits[l]->ptr) = (*(header_bits[l]->ptr) & ~1) | lsb;
      }
    }
  }

  /* Output image */
  FILE* output = fopen("./output.jpg", "w+b");

  struct jpeg_compress_struct out_cinfo;

  out_cinfo.err = jpeg_std_error(&jerr);
  out_cinfo.err->error_exit = on_error_exit;
  out_cinfo.err->output_message = on_output_message;
  jpeg_create_compress(&out_cinfo);

  jpeg_stdio_dest(&out_cinfo, output);

  out_cinfo.image_width = in_cinfo.image_width;
  out_cinfo.image_height = in_cinfo.image_height;
  out_cinfo.in_color_space = in_cinfo.jpeg_color_space;
  jpeg_set_defaults(&out_cinfo);

  jpeg_set_quality(&out_cinfo, detect_quality(&in_cinfo), TRUE);

  jpeg_write_coefficients(&out_cinfo, coeffs);
  (void)jpeg_finish_compress(&out_cinfo);

  puts("File written as output.jpg");

  bitmap_free(&bm);
  jpeg_finish_decompress(&in_cinfo);
  jpeg_destroy_decompress(&in_cinfo);
  fclose(input); 

  free(key);
  free(key_prefix);
}

