# Installing HiPhase
## Instructions
1. Navigate to the [latest release](https://github.com/PacificBiosciences/HiPhase/releases/latest) and download the tarball file (e.g. `hiphase-{version}-x86_64-unknown-linux-gnu.tar.gz`).
2. Decompress the tar file.
3. (Optional) Verify the md5 checksum.
4. Test the binary file by running it with the help option (`-h`).
5. Visit the [User guide](./user_guide.md) for details on running HiPhase.

## Example with v0.7.2
```bash
# modify this to update the version
VERSION="v0.7.2"
# get the release file
wget https://github.com/PacificBiosciences/HiPhase/releases/download/${VERSION}/hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
# decompress the file into folder ${VERSION}
tar -xzvf hiphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
cd hiphase-${VERSION}-x86_64-unknown-linux-gnu
# optional, check the md5 sum
md5sum -c hiphase.md5
# execute help instructions
./hiphase -h
```