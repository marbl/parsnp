 #!/bin/bash
make LDADD=-lMUSCLE-3.7
./bootstrap.sh
./configure
make install

