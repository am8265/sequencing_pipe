all: GoButtonAlign

GoButtonAlign: GoButtonAlign.cpp
	g++ -o GoButtonAlign GoButtonAlign.cpp -std=c++0x -I/nfs/goldstein/software/sequencing_pipe/master/sequencing_pipe/external/rapidjson/include `mysql_config --cflags --libs` -Wl,-E -Wl,-rpath,/usr/lib64/perl5/CORE -fstack-protector -L/usr/lib64/perl5/CORE -lperl -lresolv -lnsl -ldl -lm -lcrypt -lutil -lpthread -lc -D_REENTRANT -D_GNU_SOURCE -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/lib64/perl5/CORE -lcurl 