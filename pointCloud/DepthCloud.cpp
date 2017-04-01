// Jadon Wagstaff
// March 2017
/* 
   DepthCloud
   
   A grid on a globe is created with points spaced
   1/RESOLUTION size apart in the latitude direction.
   In the longitude direction the  points are spread
   1/RESOLUTION at the equator and decrease to 0 at the
   poles.
   
   For each polygon entering the stream, the id is found,
   then the corresponding evaluation is found. Then, for 
   each overlapping point on the grid, both the id and the
   evaluation is added to a vector associated with that point.
   
   Evaluation id: DD-0 LC-1 NT-2 VU-3 EN-4 CR-5
   
   Only points with depth at least one are recorded. Output 
   is a geojson file with points corresponding to 
   the grid. Secondary output is a csv file with
   latitude,longitude,id[0],status[0],id[1],status[1],...
   on each line.
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

// global constants
double RESOLUTION = 10;
double PI = 3.141592653589793;





//***** function declerations *****

// Reads in a csv with id and corresponding status
void read_status( ifstream& in, vector<int>& id, vector< vector<char> >& st );

// Markup label reader
// Returns a string with all characters until stopped by '>' character
string label_reader( ifstream& in );

// Reads in all coordinates for a polygon
// Returns a vector of (longitude, latitude) points
vector< vector<double> > coordinate_reader( ifstream& in );

// Reads in a double value from a file 
double double_reader( ifstream& in );

// Returns an integer value from a file
int int_reader( ifstream& in );

// Polygon crossings finder
// Input a polygon, and the resolution
// Return a vector of longitude values where the polygon is crossed, in descending order
vector< vector<double> >find_crossings( vector< vector< vector<double> > > polygon, double minLon, double maxLon );




int main()
{
	//***** constants and variables *****
	
	// variables
	vector< vector< vector<int> > > grid;
	vector< vector<int> > row;
	vector<int> label;
	
	vector< vector< vector<double> > > polygonWithHoles;
	vector< vector<double> > crossings;
	vector<double> coordinate;
	double minLat, maxLat, minLon, maxLon;
	int bracket = 0, count = 0;
	bool multi;
	
	char ch;
	string st;
	int id, status;
	vector<int> idLabels;
	vector< vector<char> > statusLabels;
	
	
	ifstream in;
	ofstream out;
	
	
	
	
	
	//***** create point cloud *****
	
	label.clear();
	
	// point cloud for globe grid[latitude][longitude][label]
	for( int i = -90*RESOLUTION; i < 90*RESOLUTION; i++ ) {
		row.clear();
		for( int j = floor(-180*RESOLUTION*(cos(((double(i)/RESOLUTION)*PI)/180))); j < ceil(180*RESOLUTION*(cos(((double(i)/RESOLUTION)*PI)/180))); j++ ) {
		//for( int j = -180*RESOLUTION; j < 180*RESOLUTION; j++ ) {
			row.push_back(label);
		}
		grid.push_back(row);
	}
	
	in.open("as.csv");
	read_status(in, idLabels, statusLabels);
	in.close();
	
	
	
	
	
	//***** process polygons *****
	in.open("amphibians.geojson");
	
	// loop through and find each polygon
	while( !in.eof() ) {
	
		// ensure in pointer is at correct location to read a quotation label
		while( in.peek() != '\"' && !in.eof() ) {
			in.get(ch);
		}
		
		if( !in.eof() ) {		
			st = label_reader(in);
		}
		else { st = ""; }
		
		// read in id and status
		if( st == "id_no" ) {
			in.get(ch);
			while( in.peek() == ' ' ) {
				in.get(ch);
			}
			id = int_reader(in);
			count++;
			cout << count << endl;
			for (int i = 0; i < idLabels.size(); i++) {
				if (id == idLabels[i]) {
					if (statusLabels[i][0] == 'D') {
						status = 0;
					}
					else if (statusLabels[i][0] == 'L') {
						status = 1;
					}
					else if (statusLabels[i][0] == 'N') {
						status = 2;
					}
					else if (statusLabels[i][0] == 'V') {
						status = 3;
					}
					else if (statusLabels[i][0] == 'E') {
						status = 4;
					}
					else if (statusLabels[i][0] == 'C') {
						status = 5;
					}
					else {
						cout << "Error: unidentified " << id << endl;
						status = 0;
					}
					break;
				}
				if (i == idLabels.size() - 1) {
					cout << "Error: missing " << id << endl;
					status = 0;
				}
			}
		}
		
		// read in a polygon
		else if( st == "coordinates" ) {
			in.get();
			bracket = 0;
			
			// determine whether coordinates represent Polygon or MultiPolygon
			while ( in.peek() == ' ' || in.peek() == '[' || in.peek() == ',' ) {
				in.get(ch);
				if( ch == '[' ) {
					bracket++;
				}
			}
			if( bracket == 4 ) {
				multi = true;
			}
			else {
				multi = false;
			}
			
			
			
			do {
				// find coordinates of each polygon with holes
				polygonWithHoles.clear();
				bool go = true;
				
				
				do {
					polygonWithHoles.push_back(coordinate_reader(in));
					bracket--;
					while ( in.peek() == ' ' || in.peek() == '[' || in.peek() == ']' || in.peek() == ',' ) {
						in.get(ch);
						if( ch == '[' ) {
							bracket++;
						}
						if( ch == ']' ) {
							bracket--;
						}
						if( bracket == 1 && multi == true) {
							go = false;
						}
					}
				} while ( bracket > 0  && go == true );
				
				while ( in.peek() == ' ' || in.peek() == '[' || in.peek() == ']' || in.peek() == ',' ) {
					in.get(ch);
					if( ch == '[' ) {
						bracket++;
					}
					if( ch == ']' ) {
						bracket--;
					}
				}
				
				
				
				// process each the polygon or each polygon in the multipolygon
				// find mins and maxes of polygon
				minLon = polygonWithHoles[0][0][0];
				maxLon = polygonWithHoles[0][0][0];
				minLat = polygonWithHoles[0][0][1];
				maxLat = polygonWithHoles[0][0][1];
				for( unsigned int i = 0; i < polygonWithHoles[0].size(); i++ ) {
					if( polygonWithHoles[0][i][0] < minLon ) { minLon = polygonWithHoles[0][i][0]; }
					if( polygonWithHoles[0][i][0] > maxLon ) { maxLon = polygonWithHoles[0][i][0]; }
					if( polygonWithHoles[0][i][1] < minLat ) { minLat = polygonWithHoles[0][i][1]; }
					if( polygonWithHoles[0][i][1] > maxLat ) { maxLat = polygonWithHoles[0][i][1]; }
				}
				minLat = floor( minLat * RESOLUTION ) / RESOLUTION;
				maxLat = ceil( maxLat * RESOLUTION ) / RESOLUTION;
				minLon = floor( minLon * RESOLUTION ) / RESOLUTION;
				maxLon = ceil( maxLon * RESOLUTION ) / RESOLUTION;
				
				// get the crossings in the polygon
				crossings = find_crossings( polygonWithHoles, minLat, maxLat );
				
				// find where the latitude begins in the polygon
				for( unsigned int k = 0; k < grid.size(); k++ ) {
					if( minLat <= round((-90 + k/RESOLUTION)*RESOLUTION)/RESOLUTION ) {
						
					
					
						// process all latitudes in the polygon
						for( unsigned int i = k; -90 + i/RESOLUTION  <= maxLat && i < grid.size(); i++ ) {
							// find where the longitude begins in the polygon
							for( unsigned int l = 0; l < grid[i].size(); l++ ) {
								if( minLon <= -180 + 360*(double(l)/grid[i].size()) ) {
									
								
								
									// find all crossing points within the polygon on i latitude and test them against the longitude
									for( unsigned int j = l; -180 + 360*(double(j)/grid[i].size()) <= maxLon && j < grid[i].size() && crossings.size() != 0; j++ ) {
										// find out if this point crosses the polygon
										if( crossings[i-k].size() != 0 ) {
											while( crossings[i-k].size() != 0 && crossings[i-k][0] <= -180 + 360*(double(j)/grid[i].size()) ) {
												crossings[i-k].erase( crossings[i-k].begin() );
											}
											// add the id to the grid point if there are odd number of crossings left
											if( crossings[i-k].size() % 2 == 1 ) {
												grid[i][j].push_back(id);
												grid[i][j].push_back(status);
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
				
				
				
			} while ( bracket > 0 && multi == true );
		}
	}




			
	
	
	//***** output grid *****
	
	in.close();
	out.open("points.geojson");
	out << "{ \"type\": \"FeatureCollection\",\n\t\"features\": [\n";
	for( unsigned int i = 0; i < grid.size(); i++ ) {
		for( unsigned int j = 0; j < grid[i].size(); j++ ) {
			if (grid[i][j].size() > 0) {
				out << "\t\t{ \"type\": \"Feature\",\n";
				out << "\t\t\t\"geometry\": {\"type\": \"Point\", \"coordinates\": [ ";
				out << -180 + 360*(double(j)/grid[i].size());
				out << ", ";
				out << -90 + i/RESOLUTION;
				out << " ] },\n";
				out << "\t\t\t\"properties\": {\"Species\": [ ";
				for (unsigned int k = 0; k < grid[i][j].size(); k += 2) {
					out << "[ " << grid[i][j][k] << ", " << grid[i][j][k+1] << "], ";
				}
				out << "] }\n";
				out << "\t\t},\n";
			}
		}
	}
	out << "\t]\n";
	out << "}";

	out.close();
	
	out.open("points.csv");
	out << "LATITUDE,LONGITUDE,ID,STATUS..." << endl;
	for( unsigned int i = 0; i < grid.size(); i++ ) {
		for( unsigned int j = 0; j < grid[i].size(); j++ ) {
			if (grid[i][j].size() != 0) {
				out << -90 + i/RESOLUTION << ",";
				out << -180 + 360*(double(j)/grid[i].size());
				for (unsigned int k = 0; k < grid[i][j].size(); k += 2) {
					out << "," << grid[i][j][k] << "," << grid[i][j][k+1];
				}
				out << endl;
			}
		}
	}
	
	out.close();
	
	
	
	return 0;
}















//***** functions *****

void read_status( ifstream& in, vector<int>& id, vector< vector<char> >& st ) {
	char ch;
	vector<char> chstr; int l = 0;
	
	while( !in.eof() ) {
		id.push_back(int_reader(in));
		chstr.clear();
		in.get(ch);
		while( !in.eof() && ch != '\n' ) {
			if( ch != ',' ) {
				chstr.push_back(ch);
			}
			in.get(ch);
		}
		st.push_back(chstr);
	}
}



string label_reader( ifstream& in ) {

	// variables
	char ch;
	string s;
	
	in.get(ch);
	in.get(ch);
	
	while( ch != '\"' ) {
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
	
	while( in.peek() != ']' ) {
		
		// read in longitude
		while( in.peek() == ' ' || in.peek() == '[' || in.peek() == ',' ) {
			in.get();
		}
		coordinate[0] = double_reader(in);
		
		// read in latitude
		while( in.peek() == ',' || in.peek() == ' ' ) {
			in.get();
		}
		coordinate[1] = double_reader(in);
		
		coordinates.push_back(coordinate);
		
		while( in.peek() == ' ' || in.peek() == ',' ) {
			in.get();
		}
		if( in.peek() == ']' ) {
			in.get();
		}
		while( in.peek() == ' ' || in.peek() == ',' ) {
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





vector< vector<double> >find_crossings( vector< vector< vector<double> > > polygon, double minLat, double maxLat ) {
	double beginning, ending;
	double slope, intercept;
	vector< vector<double> > crossings;
	
	// initialize the crossings vector
	for( int i = 0; i < round( abs( minLat - maxLat ) * RESOLUTION ) + 1; i++ ) {
		vector<double> row;
		crossings.push_back(row);
	}
	
	
	for( unsigned int k = 0; k < polygon.size(); k++ ) {
		for( unsigned int i = 0; i < polygon[k].size() - 1; i++ ) {
		
			// find latitude and longitude for beginning and ending
			if( polygon[k][i][1] < polygon[k][i+1][1] ) {
				beginning = ceil( polygon[k][i][1] * RESOLUTION ) / RESOLUTION;
				ending = floor( polygon[k][i+1][1] * RESOLUTION ) / RESOLUTION;
				if (beginning == polygon[k][i][1]) {
					beginning = beginning + 1/RESOLUTION;
				}
			}
			else if ( polygon[k][i][1] > polygon[k][i+1][1] ) {
				beginning = ceil( polygon[k][i+1][1] * RESOLUTION ) / RESOLUTION;
				ending = floor( polygon[k][i][1] * RESOLUTION ) / RESOLUTION;
				if (beginning == polygon[k][i+1][1]) {
					beginning = beginning + 1/RESOLUTION;
				}
			}
			
			
		
		
			if( polygon[k][i][1] != polygon[k][i+1][1] && ending >= beginning ) {
		
				// find intersections for longitude and include in crossings
				slope = ( polygon[k][i][0] - polygon[k][i+1][0] ) / ( polygon[k][i][1] - polygon[k][i+1][1] );
				intercept = polygon[k][i][0] - slope * polygon[k][i][1];
			
				for( int j = 0; j <= round(abs( beginning - ending ) * RESOLUTION); j++ ) {
					crossings[round((beginning - minLat)*RESOLUTION + j)].push_back( slope * ( beginning + j / RESOLUTION ) + intercept );
				}
			
			}
		
		}
	}
	
	for( unsigned int i = 0; i < crossings.size(); i++ ) {
		if( crossings[i].size() != 0 ) {
			sort( crossings[i].begin(), crossings[i].end() );
		}
		if( crossings[i].size() % 2 != 0 ) {
			cout << "Error " << i << endl;
		}
	}
	
	return crossings;
	
}

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
