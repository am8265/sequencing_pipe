#include <tinyxml2.h>

int main(int argc, char **argv) {

// g++ -I tinyxml2 % tinyxml2/tinyxml2.cpp && ./a.out /nfs/seqscratch1/Runs/december/180330_A00116_0093_BH7TJJDMXX/RunParameters.xml
    // using namespace 
    tinyxml2::XMLDocument doc;
    doc.LoadFile( *(argv+1) );

    // Structure of the XML file:
    // - Element "PLAY"      the root Element, which is the 
    //                       FirstChildElement of the Document
    // - - Element "TITLE"   child of the root PLAY Element
    // - - - Text            child of the TITLE Element

    // Navigate to the title, using the convenience function,
    // with a dangerous lack of error checking.
    const char* title = doc.FirstChildElement( "RunParameters" )->FirstChildElement( "RtaVersion" )->GetText();
    printf( "Name of play (1): %s\n", title );

    // Text is just another Node to TinyXML-2. The more
    // general way to get to the XMLText:
    // tinyxml2::XMLText* textNode = doc.FirstChildElement( "PLAY" )->FirstChildElement( "TITLE" )->FirstChild()->ToText();
    // title = textNode->Value();
    // printf( "Name of play (2): %s\n", title );
}
