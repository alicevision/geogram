HEXDOM_PLUGIN_PATH=/home/levy/Programming/graphite/plugins/OGF/hexdom/
rm -f hexdom/*
cp $HEXDOM_PLUGIN_PATH/algo/*.h hexdom/
cp $HEXDOM_PLUGIN_PATH/algo/*.cpp hexdom/
cd hexdom
for i in `ls *.h *.cpp`
do
    cat $i | sed -e 's|OGF/hexdom/algo|exploragram/hexdom|g' \
		 -e 's|OGF/hexdom/common|exploragram/basic|g' \
		 -e 's|hexdom_API|EXPLORAGRAM_API|g' > tempfile \
		 -e 's|std::cerr\(.*\)std::endl|GEO::Logger::out("HexDom") \1 std::endl|g'
    mv tempfile $i
done

