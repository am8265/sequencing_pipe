all: GoButton

GoButton: GoButton.cpp
	g++ -o GoButton GoButton.cpp -std=c++0x -I/nfs/goldstein/software/sequencing_pipe/master/sequencing_pipe/external/rapidjson/include -Iexternal/tinyxml2/tinyxml2/ -Iexternal/InterOp-1.1.4-Linux-GNU/include/ -I/usr/include/mysql/  -I/usr/lib64/perl5/CORE/ -Iexternal/myhtml/source/ `mysql_config --cflags --libs` -Wl,-E -Wl,-rpath,/usr/lib64/perl5/CORE -fstack-protector -L/usr/lib64/perl5/CORE -lperl -lresolv -lnsl -ldl -lm -lcrypt -lutil -lpthread -lc -D_REENTRANT -D_GNU_SOURCE -fno-strict-aliasing -pipe -fstack-protector -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/lib64/perl5/CORE -lcurl -Lexternal/myhtml/lib -lmyhtml-4.0 -Lexternal/tinyxml2/build/ -ltinyxml2 -Lexternal/InterOp-1.1.4-Linux-GNU/lib64/ -lc_csharp_run_metrics -Lexternal/InterOp-1.1.4-Linux-GNU/lib64 -lc_csharp_summary -Wl,-rpath,external/myhtml/lib -Wl,-rpath,external/tinyxml2/build/ -Wl,-rpath,external/InterOp-1.1.4-Linux-GNU/lib64/

