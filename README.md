# den2bin
Tool that converts VASP density file (CHGCAR/LOCPOT) to a binary file

## Compilation
den2bin requires the following packages
* libboost
* libbz2-dev
* libtclap-dev

Compilation is done using CMake
```
mkdir build
cd build
cmake ../src
make -j5
```

## Usage
```
./den2bin -i <DENSITY_FILENAME> -o <BINARY_FILENAME> -m <MESSAGE_HEADER>
```

Example:
```
./den2bin -i CHGCAR -o chgcar.bin -m "My special CHGCAR file"
```

## Performance
Typically, the final binary file has a compression ratio of 10-12. For example, a 300 MB file will be reduced to less than 30 MB. The compression is done in less than 5 seconds.
