# Installing HiPhase
## From conda
The easiest way to install HiPhase is through [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html):

```bash
# create a brand new conda environment and install latest HiPhase
conda create -n hiphase -c bioconda hiphase
# OR install latest into current conda environment
conda install hiphase
# OR install a specific version into current conda environment
conda install hiphase=0.8.1
```

## From GitHub
Conda updates usually lag the GitHub release by a couple days.
Use the following instructions to get the most recent version directly from GitHub:

1. Navigate to the [latest release](https://github.com/PacificBiosciences/HiPhase/releases/latest) and download the tarball file (e.g. `hiphase-{version}-x86_64-unknown-linux-gnu.tar.gz`).
2. Decompress the tar file.
3. (Optional) Verify the md5 checksum.
4. Test the binary file by running it with the help option (`-h`).
5. Visit the [User guide](./user_guide.md) for details on running HiPhase.

### Example with v0.7.2
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