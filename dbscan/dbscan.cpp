// Jadon Wagstaff
// March 2017
/* dbscan algorithm is implemented to cluster
   centroid data of polygons based on density
   */

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>
#include<algorithm>
#include<vector>
#include<cmath>

using namespace std;

// global variables
vector< vector<double> > points;
vector< vector< vector<int> > > distances;

// calculates haversine distance between two (longitude, latitude) points in kilometers
double Haversine( vector<double> p1, vector<double> p2 );

// subroutine of dbscan, used to choose new labels
void NbhdJump( int p, int label, int d, int minPts );

// Comma sepparated value reader
// Returns a string with all characters until stopped by ',' character
void CsvReader( ifstream& in );

// Reads in a double value from a file 
double DoubleReader( ifstream& in );


int main()
{
	//***** constants and variables *****
	
	// variables
	
	vector<int> d;
	int cluster = 0;
	double epsilon = 50;
	int choosed, minPts, size;
	
	d.push_back(25);
	d.push_back(50);
	d.push_back(75);
	d.push_back(100);
	d.push_back(125);
	d.push_back(150);
	
	
	
	//***** process points *****
	vector<double> point, row;
	char ch, go = 'y';
	ifstream in;
	ofstream out;
	
	// each point has a longitude, latitude, and label
	point.push_back(0);
	point.push_back(0);
	point.push_back(0);
	
	in.open("Amphibians.txt");
	
	in.get(ch);
	while( ch != '\n' ) {
		in.get(ch);
	}
	
	while( !in.eof() ) {
		for( int i = 0; i < 26; i++ ) {
			if( in.eof() ) {break;}
			CsvReader(in);
		}
		if( in.eof() ) {break;}
		point[0] = DoubleReader(in);
		point[1] = DoubleReader(in);
		points.push_back(point);
	}
	
	in.close();
	
	in.open("Reptiles.txt");
	
	in.get(ch);
	while( ch != '\n' ) {
		in.get(ch);
	}
	
	while( !in.eof() ) {
		for( int i = 0; i < 26; i++ ) {
			if( in.eof() ) {break;}
			CsvReader(in);
		}
		if( in.eof() ) {break;}
		point[0] = DoubleReader(in);
		point[1] = DoubleReader(in);
		points.push_back(point);
	}
	
	in.close();
	
	in.open("TerrestrialMammals.txt");
	
	in.get(ch);
	while( ch != '\n' ) {
		in.get(ch);
	}
	
	while( !in.eof() ) {
		for( int i = 0; i < 26; i++ ) {
			if( in.eof() ) {break;}
			CsvReader(in);
		}
		if( in.eof() ) {break;}
		point[0] = DoubleReader(in);
		point[1] = DoubleReader(in);
		points.push_back(point);
	}
	
	in.close();
	
	random_shuffle( points.begin(), points.end() );
	
	cout << "Points: " << points.size() << endl;
	
	
	
	//***** process distances *****
	// create an array of linked lists with the indeces of each point close to a point
	vector<int> blankRow;
	vector< vector<int> > blankD;
	double havD;
	
	for( int i = 0; i < points.size(); i++ ) {
		blankD.push_back(blankRow);
	}
	for( int k = 0; k < d.size(); k++ ) {
		distances.push_back(blankD);
	}
	
	for( int i = 0; i < points.size() - 1; i++ ) {
		for( int j = i + 1; j < points.size(); j++ ) {
			havD = Haversine( points[i], points[j] );
			if( havD < d[0] ) {
				distances[0][i].push_back(j);
				distances[0][j].push_back(i);
			}
			else {
				for( int k = 1; k < d.size(); k++ ) {
					if( havD >= d[k-1] && havD < d[k] ) {
						distances[k][i].push_back(j);
						distances[k][j].push_back(i);
					}
				}
			}
		}
		if( i - floor( double(i)/1000 )*1000 == 1 ) {cout << i << endl;}
	}
	
	
	
	
	
	
	//***** dbscan algorithm *****
	while( go != 'n' ) {
		
		// this allows you to choose new epsilon and minPts values without recalculating all distances
		do {
			cout << "D?";
			cin >> choosed;
		} while ( choosed > 6 );
		cout << "Points?";
		cin >> minPts;
		
		// initalize labels
		cluster = 0;
		for( unsigned int i = 0; i < points.size(); i++ ) {
			points[i][2] = 0;
		}
		
		// dbscan
		for( unsigned int i = 0; i < points.size(); i++ ) {
			// find unlabeled points
			if( points[i][2] == 0 ) {
				// find number of points within epsilon of a point
				size = 0;
				for( int k = 0; k < choosed; k++ ) {
					size += distances[k][i].size();
				}
				// if it is a center, create a cluster around it
				if( size >= minPts ) {
					cluster++;
					points[i][2] = cluster;
					NbhdJump( i, cluster, choosed, minPts );
				}
			}
		}
		
		// output results
		cout << cluster << endl;
		
		out.open("results.geojson");	
		out << "{ \"type\": \"FeatureCollection\",\n\t\"features\": [\n";
		for( int k = 1; k <= cluster; k++ ) {
			out << "\t\t{ \"type\": \"Feature\",\n";
			out << "\t\t\t\"geometry\": {\"type\": \"MultiPoint\", \"coordinates\": [ ";
			bool first = true;
			for( int i = 0; i < points.size(); i++ ) {
				if( points[i][2] == k ) {
					if( first != true ) {out << ", ";}
					first = false;
					out << "[" << points[i][0] << ", " << points[i][1] << "]";
				}
			}
			out << " ] },\n";
			out << "\t\t\t\"properties\": {\"cluster\": " << k << "}\n";
			out << "\t\t}";
			if( k != cluster ) {out << ",";}
			out << "\n";
		}
		out << "\t]\n";
		out << "}";
	
		out.close();
		
		// loop exit strategy
		cout << "Continue?(y/n)";
		cin >> go;
	}
	
	
	
	return 0;
}

void NbhdJump( int p, int label, int d, int minPts ) {
	int size = 0;
	for( int k = 0; k < d; k++ ) {
		for( int i = 0; i < distances[k][p].size(); i++ ) {
			if( points[distances[k][p][i]][2] != label) {
				points[distances[k][p][i]][2] = label;
				// if this point has large enough size jump to its neighborhood
				size = 0;
				for( int l = 0; l < d; l++ ) {
					size += distances[l][distances[k][p][i]].size();
				}
				if( size >= minPts ) {
					NbhdJump( distances[k][p][i], label, d, minPts );
				}
			}
		}
	}
}




double Haversine( vector<double> p1, vector<double> p2 ) {
	const double PI = 3.141592653589793;
	double lon1, lon2, lat1, lat2;
	lon1 = (p1[0]*PI)/180;
	lon2 = (p2[0]*PI)/180;
	lat1 = (p1[1]*PI)/180;
	lat2 = (p2[1]*PI)/180;
	
	
	const double R = 6366;
	double a, c, d;
	a = pow( sin((lat1 - lat2)/2), 2 ) + cos(lat1)*cos(lat2)*pow( sin((lon1 - lon2)/2), 2 );
	c = 2 * atan2(sqrt(a), sqrt(1-a));
	d = R * c;
	return d;
}





void CsvReader( ifstream& in ) {
	char ch;
	
	in.get(ch);
	
	while( ch != ',' && !in.eof() ) {
		if( ch == '"' ) {
			in.get(ch);
			while( ch != '"' && !in.eof() ) {
				in.get(ch);
			}
		}
		in.get(ch);
	}
	
}








double DoubleReader( ifstream& in ) {
	char ch;
	vector<char> adec, bdec, expc;
	double val = 0, exp = 0;
	bool neg, small;	
	
	
	// read in value
	if (in.peek() == '-')
	{
		neg = true;
		in.get(ch);
	}
	else {neg = false;}
	
	// find value before decimal
	while( in.peek() == '0' || in.peek() == '1' || in.peek() == '2' || in.peek() == '3' || in.peek() == '4' || in.peek() == '5' || in.peek() == '6' || in.peek() == '7' || in.peek() == '8' || in.peek() == '9' ) {
		in.get(ch);
		bdec.push_back(ch);
	}
	
	// find value after decimal
	if ( in.peek() == '.' ) {
		in.get(ch);
		while ( in.peek() == '0' || in.peek() == '1' || in.peek() == '2' || in.peek() == '3' || in.peek() == '4' || in.peek() == '5' || in.peek() == '6' || in.peek() == '7' || in.peek() == '8' || in.peek() == '9' ) {
			in.get(ch);
			adec.push_back(ch);
		}
	}
	
	// find exponential value
	if ( in.peek() == 'e' ) {
		in.get(ch);
		if ( in.peek() == '-' ) {
			small = true;
			in.get(ch);
		}
		else { small = false; }
		while ( in.peek() == '0' || in.peek() == '1' || in.peek() == '2' || in.peek() == '3' || in.peek() == '4' || in.peek() == '5' || in.peek() == '6' || in.peek() == '7' || in.peek() == '8' || in.peek() == '9' ) {
			in.get(ch);
			expc.push_back(ch);
		}
	}
	

	// add value before decimal to double
	for ( unsigned int i = 0; i < bdec.size(); i++ ) {
		val = val + double((int)(bdec[i] - '0'))*(pow(double(10), (bdec.size() - i - 1)));
	}
	
	// add value after decimal to double
	if (adec.empty() == false)
	{
		for (unsigned int i = 0; i < adec.size(); i++)
		{
			val = val + double((int)(adec[i] - '0'))*(pow(double(10), (-1*(double(i) + 1))));
		}
	}
	
	// add exponential value to double
	if (expc.empty() == false)
	{
		for (unsigned int i = 0; i < expc.size(); i++)
		{
			exp = exp + (int)(expc[i] - '0')*(pow(10, (expc.size() - i - 1)));
		}
		if (small == true) 
		{
			exp = exp*int(-1);
		}
		val = val*(pow(10, exp));
	}
	
	// make double negative if necessary
	if (neg == true)
	{
		val = val*double(-1);
	}
	
	in.get(ch);
	
	
	return val;
						
}































































