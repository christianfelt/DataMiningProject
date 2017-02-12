#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>
#include<algorithm>
#include<vector>
#include<cmath>

using namespace std;


// Utah coordinates
// North:  42, -114     42, -111
// Center: 41, -111     41, -109    
// South:  37, -109     37, -114

// Rectangular overview of Utah coordinates
// nw: 45, -117 ne: 45, -104
// sw: 36, -117 se: 36, -104

// global constants
double N = 45, S = 36, E = -104, W = -117;
double RESOLUTION = 5;

// Markup label reader
// Returns a string with all characters until stopped by '>' character
string ml_reader( ifstream& in );

// Reads in all coordinates for a polygon
// Returns a vector of (latitude, longitude) points
vector< vector<double> > coordinate_reader( ifstream& in );

// Reads in a double value from a file 
double double_reader( ifstream& in );

// Returns an integer value from a file
int int_reader( ifstream& in );

// Polygon crossings finder
// Input a polygon, and the resolution
// Return a vector of longitude values where the polygon is crossed, in descending order
vector< vector<double> > find_crossings( vector< vector<double> > polygon );




int main()
{
	//***** constants and variables *****
	
	// variables
	vector< vector< vector<int> > > grid;
	vector< vector<int> > row;
	vector<int> label;
	
	vector< vector<double> > polygon, crossings;
	vector<double> coordinate;
	
	char ch;
	string st;
	int id;
	double minLat, maxLat, minLon, maxLon;
	ifstream in;
	
	
	//***** create point cloud *****
	
	label.clear();
	
	// point cloud for rectangle
	for( int i = 0; i <= (N - S) * RESOLUTION; i++ ) {
		for( int j = 0; j <= (E - W) * RESOLUTION; j++ ) {
			row.push_back(label);
		}
		grid.push_back(row);
		row.clear();
	}
	
	
	//***** process polygons *****
	
	in.open("test.kml");
	
	
	// loop through and find each polygon
	while( !in.eof() ) {
	
		// ensure in pointer is at correct location to read a new markup label
		while( in.peek() != '<' && !in.eof() ) {
			in.get(ch);
		}
		
		if( !in.eof() ) {		
			st = ml_reader(in);
		}
		else { st = ""; }
		
		// read in label and make a new id or read in a polygon
		if( st == "SimpleData name=\"id_no\"" ) {
			while( in.peek() == ' ' ) {
				in.get(ch);
			}
			id = int_reader(in);
		}
		else if( st == "coordinates" ) {
			
			polygon = coordinate_reader(in);
			
			// find mins and maxes of polygon
			minLat = polygon[0][0];
			maxLat = polygon[0][0];
			minLon = polygon[0][1];
			maxLon = polygon[0][1];
			for( unsigned int i = 0; i < polygon.size(); i++ ) {
				if( polygon[i][0] < minLat ) { minLat = polygon[i][0]; }
				if( polygon[i][0] > maxLat ) { maxLat = polygon[i][0]; }
				if( polygon[i][1] < minLon ) { minLon = polygon[i][1]; }
				if( polygon[i][1] > maxLon ) { maxLon = polygon[i][1]; }
			}
			if( minLat < S && maxLat > S ) {
				minLat = S;
			}
			if( maxLat > N && minLat < N ) {
				maxLat = N;
			}
			minLat = ceil( minLat * RESOLUTION ) / RESOLUTION;
			maxLat = floor( maxLat * RESOLUTION ) / RESOLUTION;
			minLon = ceil( minLon * RESOLUTION ) / RESOLUTION;
			maxLon = floor( maxLon * RESOLUTION ) / RESOLUTION;
			
			// compare polygons to grid
			// check to see if polygon overlaps search area
			if( ( ( S < minLat && minLat < N )  || ( S < maxLat && maxLat < N ) ) && ( ( W < minLon && minLon < E ) || ( W < maxLon && maxLon < E ) ) ) {
				// get the crossings in the polygon
				crossings = find_crossings( polygon );
			
				// find where the latitude begins in the polygon
				for( unsigned int k = 0; k < grid.size(); k++ ) {
					if( minLat <= S + k / RESOLUTION ) {
						
					
					
						// process all latitudes in the polygon
						for( unsigned int i = k; S + i / RESOLUTION <= maxLat && i < grid.size(); i++ ) {
							// find where the longitude begins in the polygon
							for( unsigned int l = 0; l < grid[i].size(); l++ ) {
								if( minLon <= W + l / RESOLUTION  ) {
								
								
								
									// find all crossing points within the polygon on i latitude and test them against the longitude
									for( unsigned int j = l; W + j / RESOLUTION <= maxLon && j < grid[i].size() && crossings.size() != 0; j++ ) {
										// find out if this point crosses the polygon
										if( crossings[i-k].size() != 0 ) {
											if( crossings[i-k][0] <= W + j / RESOLUTION ) {
												crossings[i-k].erase( crossings[i-k].begin() );
											}
											// add the id to the grid point if there are odd number of crossings left
											if( crossings[i-k].size() % 2 == 1 ) {
												grid[i][j].push_back(id);
											}
										}
									} // longitude
									
									
									
									break;
								}
							}
						} // latitude
						
						
						
						break;
					}
				}
										
			} // checking polygon against grid	
			
			
			
			
		} // processing polygon
		
		
	} // looping through polygons
			
	
	
	// TODO output grid
	
	for( unsigned int i = 0; i < grid.size(); i++ ) {
		for( unsigned int j = 0; j < grid[i].size(); j++ ) {
			if( grid[i][j].size() != 0 ) {
				for( unsigned int k = 0; k < grid[i][j].size(); k++ ) {
					cout  << " " << grid[i][j][k];
				}
			}
			cout << ",";
		}
		cout << endl;
	}
	
	
	
	return 0;
}















string ml_reader( ifstream& in ) {

	// variables
	char ch;
	string s;
	
	in.get(ch);
	in.get(ch);
	
	while( ch != '>' ) {
		s.push_back(ch);
		in.get(ch);
	}
	
	return s;

}


vector< vector<double> > coordinate_reader( ifstream& in ) {
	vector< vector<double> > coordinates;
	vector<double> coordinate;
	
	coordinate.push_back(0);
	coordinate.push_back(0);
	
	while( in.peek() != '<' ) {
		
		// read in latitude
		while( in.peek() == ' ' ) {
			in.get();
		}
		if( in.peek() != '<' ) {
			coordinate[0] = double_reader(in);
		}
		
		// read in latitude
		while( in.peek() == ',' || in.peek() == ' ' ) {
			in.get();
		}
		if( in.peek() != '<' ) {
			coordinate[1] = double_reader(in);
		}
		
		coordinates.push_back(coordinate);
		
		while( in.peek() == ' ' ) {
			in.get();
		}
		
	}
	
	return coordinates;
	
}
	




double double_reader( ifstream& in ) {
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
	
	
	return val;
						
}

int int_reader( ifstream& in ) {

	// variables
	char ch;
	vector<char> chstr;
	int num = 0;
	
	// get number as string
	while( in.peek() == '0' || in.peek() == '1' || in.peek() == '2' || in.peek() == '3' || in.peek() == '4' || in.peek() == '5' || in.peek() == '6' || in.peek() == '7' || in.peek() == '8' || in.peek() == '9' ) {
		in.get(ch);
		chstr.push_back(ch);
	}
	
	// convert string to integer
	for (unsigned int i = 0; i < chstr.size(); i++)
	{
		num = num + (int)(chstr[i] - '0')*(pow(double(10), (chstr.size() - i - 1)));
	}
	
	return num;

}


vector< vector<double> >find_crossings( vector< vector<double> > polygon ) {
	double minLat, maxLat;
	double beginning, ending;
	double slope, intercept;
	vector< vector<double> > crossings;
	
	// find mins and maxes of polygon
	minLat = polygon[0][0];
	maxLat = polygon[0][0];
	for( unsigned int i = 0; i < polygon.size(); i++ ) {
		if( polygon[i][0] < minLat ) { minLat = polygon[i][0]; }
		if( polygon[i][0] > maxLat ) { maxLat = polygon[i][0]; }
	}
	if( minLat < S && maxLat > S ) {
		minLat = S;
	}
	if( maxLat > N && minLat < N ) {
		maxLat = N;
	}
	minLat = ceil( minLat * RESOLUTION ) / RESOLUTION;
	maxLat = floor( maxLat * RESOLUTION ) / RESOLUTION;
	
	// initialize the crossings vector
	for( int i = 0; i < floor( abs( minLat - maxLat ) * RESOLUTION ) + 1; i++ ) {
		vector<double> row;
		crossings.push_back(row);
	}
	
	for( unsigned int i = 0; i < polygon.size() - 1; i++ ) {
		
		// find latitude and longitude for beginning and ending
		if( polygon[i][0] < polygon[i+1][0] ) {
			beginning = ceil( polygon[i][0] * RESOLUTION ) / RESOLUTION;
			ending = floor( polygon[i+1][0] * RESOLUTION ) / RESOLUTION;
		}
		else if ( polygon[i][0] > polygon[i+1][0] ) {
			beginning = ceil( polygon[i+1][0] * RESOLUTION ) / RESOLUTION;
			ending = floor( polygon[i][0] * RESOLUTION ) / RESOLUTION;
		}
		
		// adjust beginning and ending if necessary
		if( beginning < minLat && ending > minLat ) {
			beginning = minLat;
		}
		if( ending > maxLat && beginning < maxLat ) {
			ending = maxLat;
		}
		
		if( beginning <= ending && beginning >= minLat && ending <= maxLat && polygon[i][0] != polygon[i+1][0] ) {
		
			// find intersections for longitude and include in crossings
			slope = ( polygon[i][1] - polygon[i+1][1] ) / ( polygon[i][0] - polygon[i+1][0] );
			intercept = polygon[i][1] - slope * polygon[i][0];
			
			for( int j = 0; j <= abs( beginning - ending ) * RESOLUTION; j++ ) {
				crossings[ abs( minLat - beginning ) * RESOLUTION + j ].push_back( slope * ( beginning + j / RESOLUTION ) + intercept );
			}
			
		}
		
	}
	
	if( crossings[0].size() % 2 != 0 ) {
		crossings[0].clear();
	}
	
	for( unsigned int i = 1; i < crossings.size(); i++ ) {
		if( crossings[i].size() != 0 ) {
			sort( crossings[i].begin(), crossings[i].end() );
		}
		if( crossings[i].size() % 2 != 0 ) {
			crossings[i] = crossings[i-1];
		}
	}
	
	return crossings;
	
}

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
