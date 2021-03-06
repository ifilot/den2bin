# den2bin
Tool that converts VASP density file (CHGCAR/LOCPOT) to a binary file

## Purpose
den2bin can compress CHGCAR/LOCPOT files either lossless or lossy. In the latter procedure, a 3D discrete cosine transformation is used wherein the higher frequency cosine coefficients are dropped (similar in the JPEG format).

## Obtaining den2bin from the repository

Install the GPG key

```
wget -O - http://www.zuidcoach.nl/repos/apt/conf/ivofilotsoftware.gpg.key|sudo apt-key add -
```
Add repository

```
echo "deb http://zuidcoach.nl/repos/apt/debian/ ivofilot-software main" | sudo tee -a /etc/apt/sources.list.d/ivofilot-software.list
```

Update apt sources and install `den2bin`

```
sudo apt-get update
sudo apt-get install den2bin
```

## Compilation
den2bin requires the following packages (tested on Debian)
* libboost
* libbz2-dev
* libglm-dev
* libtclap-dev

Compilation is done using CMake
```
mkdir build
cd build
cmake ../src
make -j5
```

## Usage

### Lossy compression
```
# to pack (note the "-l")
./den2bin -i <DENSITY_FILENAME> -o <BINARY_FILENAME> -m <MESSAGE_HEADER> -l -b <BLOCKSIZE> -q <QUALITY>

# to unpack (same as for the lossless version)
./den2bin -x -i <BINARY_FILENAME> -o <DENSITY_FILENAME>
```

Example:
```
# to pack (note the "-l")
./den2bin -i CHGCAR -o chgcar.dct -m "My special CHGCAR file" -l -b 4 -q 2

# to unpack (same as for the lossless version)
./den2bin -i chgcar.dct -o CHGCAR -x
```

### Lossless compression
```
# to pack
./den2bin -i <DENSITY_FILENAME> -o <BINARY_FILENAME> -m <MESSAGE_HEADER>

# to unpack
./den2bin -x -i <BINARY_FILENAME> -o <DENSITY_FILENAME>
```

Example:
```
# to pack
./den2bin -i CHGCAR -o chgcar.bin -m "My special CHGCAR file"

# to unpack
./den2bin -i chgcar.bin -o CHGCAR -x
```

## Performance
For the lossless compression, the final binary file has a compression ratio of 10-12. For example, a 300 MB file will be reduced to less than 30 MB. The compression is done in less than 5 seconds.

For lossy compression, the compression ratio can be as large as 250. For example, a 500 MB file will be reduced to about 2 MB. Compression times are typically between 1-10 seconds (tested on a i7-4790K using 8 threads).

## Showcase

All images have been generated using [EDP](https://github.com/ifilot/edp).

### CHGCAR files

Recommended setting for CHGCAR is `Blocksize 4, Quality 4`

![CHGCAR compression as a function of blocksize and quality](https://raw.githubusercontent.com/ifilot/den2bin/master/examples/chgcar.jpg)

### LOCPOT files

Recommended setting for LOCPOT is `Blocksize 4, Quality 6`

![LOCPOT compression as a function of blocksize and quality](https://raw.githubusercontent.com/ifilot/den2bin/master/examples/locpot.jpg)
