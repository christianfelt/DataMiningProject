// Jadon Wagstaff
// March 2017
/* This code helps fill in information missing from out data.
   It reads in our data and finds an assessment for each animal id.
   I did not include the token information for the IUCN api since it is
   not supposed to be public.
   I used partial code retrieved from http://www.cplusplus.com/forum/unices/45878/
   on march 27.
   */


#include <iostream>
#include<fstream>
#include <string>

#include <curl/curl.h>

using namespace std;

string data; //will hold the url's contents

// Comma sepparated value reader
// Returns a string with all characters until stopped by ',' character
string CsvReader( ifstream& in );

//callback must have this declaration
size_t writeCallback(char* buf, size_t size, size_t nmemb, void* up);


int main() {
    CURL* curl;
    string begin = "http://apiv3.iucnredlist.org/api/v3/species/id/";
    string end = "?token=c40f354d561a77d8110d7fc4fab248a0092e4b4f361cd0b1ac16b5cc64cbd669";
    string id, url, assessment;
    int count = 0;
    char ch;
    ifstream in;
    ofstream out;
    
    curl_global_init(CURL_GLOBAL_ALL);
    curl = curl_easy_init();
    
    in.open("data/centroids/TerrestrialMammals.txt");
    out.open("TerrestrialMammals_Status.csv");
	
	in.get(ch);
	while( ch != '\n' ) {
		in.get(ch);
	}
	
	
	while( !in.eof() ) {
		// find next id to process
		while( ch != ',' ) {
			in.get(ch);
		}
		id = CsvReader(in);
		in.get(ch);
		while( ch != '\n' && !in.eof() ) {
			in.get(ch);
		}
		url = begin + id + end;
		
		// get data from url associated with id
		data.clear();

		curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &writeCallback);
		curl_easy_perform(curl);
		
		// output results
		count++;
		assessment = data.substr( data.find("category\":\"") + 11, data.find("\",\"cri", data.find("category\":\"") + 11) - data.find("category\":\"") - 11 );
		
		cout << count << " " << id << " " << assessment  << endl;
		out << id << "," << assessment << endl;
	}
	
	curl_easy_cleanup(curl);
	curl_global_cleanup();
	in.close();
	out.close();

    

    return 0;
}


size_t writeCallback(char* buf, size_t size, size_t nmemb, void* up) { 
    for (unsigned int c = 0; c<size*nmemb; c++)
    {
        data.push_back(buf[c]);
    }
    return size*nmemb; //tell curl how many bytes we handled
}

string CsvReader( ifstream& in ) {
	char ch;
	string s;
	
	in.get(ch);
	
	while( ch != '.' && ch != ',' && !in.eof() ) {
		s.push_back(ch);
		in.get(ch);
	}
	
	return s;
	
}

