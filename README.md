This project is simple POC showing that an image could theorically have more than 2 messages hidden which can be extracted using Outguess.

It bruteforces a prefix in the form of `<number>_` such that the resulting image will have `n` 128-byte messages embedded in it. The keys for the messages are `<number>_0` through `<number>_<n - 1>`. This POC does not bruteforce for the message contents, so the `n` messages will be junk data.

Usage: `./n_outguess <image.jpg> <n>`

Once finished, an image named `output.jpg` will be written.

## Building
Run `make`

## Demo
The image `poc.jpg` contains 9 Outguessable messages. The keys for the `n` messages are: `165445_0` through `165445_8`.

[Demo image](https://commons.wikimedia.org/wiki/File:Honeywell_Hall_Effect_data_function_key_front.jpg), dork vader. Public Domain.

## License
This tool uses code from the Independent JPEG Group's JPEG software, see [jpeg-6b-steg/README](jpeg-6b-steg/README)

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)
