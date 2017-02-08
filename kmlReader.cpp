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

// Markup label reader
// Returns a string with all characters until stopped by '>' character
string ml_reader( ifstream& in );

// Markup content reader
// Returns a string with all characters until stopped by '<' character
string mc_reader( ifstream& in );

// TODO
// Polygon crossings finder
// Input a polygon, and a latitude
// Return a vector of longitude values where the polygon is crossed, in descending order
find_crossings( polygon, latitude );

int main()
{
	//***** constants and variables *****
	
	// constants
	double N = 45, S = 36, E = -104, W -117;
	double RESOLUTION = 1000;
	
	// variables
	vector< vector< vector<int> > > grid;
	vector< vector<int> > row;
	vector<int> label;
	
	vector< vector<double> > polygon;
	vector<double> coordinate;
	
	char ch;
	int id;
	double minLat, maxLat, minLon, maxLon;
	ifstream in;
	
	
	//***** create point cloud *****
	
	label.clear();
	
	// point cloud for rectangle
	for( int i = 0; i < N - S; i++ ) {
		for( int ir = 0; ir < RESOLUTION; ir++ ) {
			for( int j = 0; j < E - W; j++ ) {
				for( int jr = 0; jr < RESOLUTION; jr++ ) {
					row.push_back(label);
				}
			}
			grid.push_back(row);
			row.clear();
		}
	}
	
	
	//***** process polygons *****
	
	in.open(amphibians.kml);
	
	
	// loop through and find each polygon
	while( !in.eof() ) {
	
		// ensure in pointer is at correct location to read a new markup label
		do {
			in.get(ch)
		} while( ch != '<' )
		
		// read in label and make a new id or read in a polygon
		if( ml_reader(in) == "SimpleData name=\"id_no\"" ) {
			// TODO mc_reader output to number for id
		}
		else if( ml_reader(in) == "Polygon" ) {
		
			// TODO get coordinates for polygon
			
			// find mins and maxes of polygon
			minLat = polygon[0][0];
			maxLat = polygon[0][0];
			minLon = polygon[0][1];
			maxLon = polygon[0][1];
			for( unsigned int i = 0; i < polygon.size(); i++ ) {
				if( polygon[i][0] < minLat ) { minLat = polygon[i][0]; }
				if( polygon[i][0] > maxLat ) { maxLat = polygon[i][0]; }
				if( polygon[i][1] < minLon ) { minLon = polygon[i][1]; }
				if( polygon[i][1] < maxLon ) { maxLon = polygon[i][1]; }
			}
			
			
			
			// compare polygons to grid
			// check to see if polygon overlaps search area
			if( ( ( S < minLat && minLat < N )  || ( S < maxLat && maxLat < N ) ) && ( ( W < minLon && minLon < E ) || ( W < maxLon && maxLon < E ) ) ) {
			
				// find where the latitude begins in the polygon
				for( unsigned int k = 0; k < grid.size(); k++ ) {
					if( minLat <= S + RESOLUTION * k ) {
					
					
					
						// process all latitudes in the polygon
						for( unsigned int i = k; S + RESOLUTION * i > maxLon || i < grid.size(); i++ ) {
							// find where the longitude begins in the polygon
							for( unsigned int l = 0; l < grid[i].size(); l++ ) {
								if( minLon <= W + RESOLUTION * l ) {
								
								
								
									// find all crossing points within the polygon on i latitude and test them against the longitude
									vector<double> crossings = find_crossings( polygon, S + RESOLUTION + i );
									for( unsigned int j = l; W + RESOLUTION * j < maxLat || j < grid[i].size() || crossings.size() == 0; j++ ) {
										// find out if this point crosses the polygon
										if( crossings[crossings.size()-1] > W + RESOLUTION * j ) {
											crossings.pop_back();
										}
										// add the id to the grid point if there are odd number of crossings left
										if( crossings % 2 == 1 ) {
											grid[i][j].push_back(id);
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
	
	
	
	return 0;
}















string ml_reader( ifstream& in ) {

	// variables
	char ch = in.get();
	string s;
	
	if( ch == '/' ) {
		in.get(ch);
	}
	
	while( ch != '>' ) {
		s.push_back(ch);
		in.get(ch);
	}
	
	return s;

}

string mc_reader( ifstream& in ) {

	// variables
	char ch = in.get();
	string s;
	
	while( ch != '<' ) {
		s.push_back(ch);
		in.get(ch);
	}
	
	return s;

}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
