# WASM Version

This folder contains code for compiling several of the functions in the
`cmstatrExt` package to Web Assembly (WASM). This allows the functions
to be run in a web browser without requiring R or for the user to install
any software. This is used in the website that accompanies this package.

This code relies on various functions provided by R. However, R isn't used
directly. Instead, a copy of the relevant parts of the R source code are
included. The relevant parts of the R source code are in the `nmath`
folder and are adapted as necessary to work for this application.

## License
The code in this folder has the same licenses as the rest of the package.
The entire package is licensed under the AGPL-3 license. This means that
you can use this code or the resulting WASM and JS files, subject to the
AGPL-3 license.

## Building
To build, run `make`. If you wish to test an `html` file that executes
the WASM code, you'll need to run a web server from your computer
(due to browser security restrictions). You can run `make devserver`
to do so (this requires that you have `python` in your path).


## Development Setup
### Installing `emsdk`

```bash
# Get the emsdk repo
git clone https://github.com/emscripten-core/emsdk.git

# Enter that directory
cd emsdk

# Fetch the latest version of the emsdk (not needed the first time you clone)
git pull

# Download and install the latest SDK tools.
./emsdk install latest

# Make the "latest" SDK "active" for the current user. (writes .emscripten file)
./emsdk activate latest
```
