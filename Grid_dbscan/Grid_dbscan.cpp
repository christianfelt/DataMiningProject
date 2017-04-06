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

// global variables and constants
vector< vector< vector<double> > > points;
vector< vector< vector<int> > > distances;
double epsilon = 50;
int minPts = 100;
const int RESOLUTION = 5;

// calculates haversine distance between two (longitude, latitude) points in kilometers
double Haversine( vector<double> p1, vector<double> p2 );

// subroutine of dbscan, used to choose new labels
void NbhdJump( int i, int j, int cluster );

// Reads in a double value from a file 
double DoubleReader( ifstream& in );


int main()
{
	//***** constants and variables *****
	int cluster = 0;
	
		vector<double> coord1, coord2;
		coord1.push_back(0);
		coord1.push_back(0);
		coord2.push_back(double(1)/RESOLUTION);
		coord2.push_back(0);
	int ptRadius = ceil(epsilon/Haversine(coord1, coord2));
	
	
	
	//***** process points *****
	vector< vector<double> > row2;
	vector<double> point, row;
	char ch, go = 'y';
	ifstream in;
	ofstream out;
	
	// each point has a longitude, latitude, and label
	
	in.open("mammals.csv");
	cout << "Processing input..." << endl;
	
	in.get(ch);
	while( ch != '\n' ) {
		in.get(ch);
	}
	
	while( !in.eof() ) {
		row.clear();
		
		// latitude
		row.push_back(DoubleReader(in));
		in.get(ch);
		
		// latitude
		row.push_back(DoubleReader(in));
		in.get(ch);
		
		// label
		row.push_back(0);
		
		// identification numbers and coorelating statuse
		while( ch != '\n' && !in.eof() ) {
			row.push_back(DoubleReader(in));
			in.get(ch);
		}
		
		// create new vector for each new value of latitude
		if( row2.size() != 0 && row2[row2.size()-1][0] != row[0] ) {
			points.push_back(row2);
			row2.clear();
		}
		row2.push_back(row);
	}
	
	in.close();
	cout << "Input processed." << endl;
	
	
	
	//***** process distances *****
	
	
	// create an array of linked lists with the indeces of each point close to a point
	vector<int> blankRow;
	vector< vector<int> > blankRow2;
	double havD;
	int size, count;
	
	for( unsigned int i = 0; i < points.size(); i++ ) {
		for( unsigned int j = 0; j < points[i].size(); j++ ) {
			blankRow2.push_back(blankRow);
		}
		distances.push_back(blankRow2);
		blankRow2.clear();
	}
	
	for( int i = 0; i < points.size(); i++ ) {
		if( i%50 == 0 ) {
			cout << i << endl;
		}
		for( int j = 0; j < points[i].size(); j++ ) {
			// find number of points within epsilon of a point
			size = (points[i][j].size() - 3)/2;
			count = 1;
			
			// find points left of point
			for( int m = 1; m <= ptRadius && j - m >= 0; m++ ) {
				if( Haversine(points[i][j], points[i][j-m]) <= epsilon ) {
					distances[i][j].push_back(i);
					distances[i][j].push_back(j-m);
					size += (points[i][j-m].size() - 3)/2;
					count++;
				}
			}
			
			// find points right of given latitude
			for( int m = 1; m <= ptRadius && points[i].size() - j - m > 0; m++ ) {
				if( Haversine(points[i][j], points[i][j+m]) <= epsilon ) {
					distances[i][j].push_back(i);
					distances[i][j].push_back(j+m);
					size += (points[i][j+m].size() - 3)/2;
					count++;
				}
			}
			
			// find points above latitude
			for( int k = 1; k <= ptRadius && points.size() - i - k > 0; k++ ) {
				// find where points are for a given latitude
				for( int l = 1; l < points[i+k].size(); l++ ) {
					if( points[i+k][l][1] > points[i][j][1] || points[i+k].size() - 1 == l ) {
					
						// find points left of point on a given latitude
						for( int m = 0; m <= ptRadius + 1 && l - m >= 0; m++ ) {
							if( Haversine(points[i][j], points[i+k][l-m]) <= epsilon ) {
								distances[i][j].push_back(i+k);
								distances[i][j].push_back(l-m);
								size += (points[i+k][l-m].size() - 3)/2;
								count++;
							}
						}
						
						// find points right of point on a given latitude
						for( int m = 1; m <= ptRadius && points[i+k].size() - l - m > 0; m++ ) {
							if( Haversine(points[i][j], points[i+k][l+m]) <= epsilon ) {
								distances[i][j].push_back(i+k);
								distances[i][j].push_back(l+m);
								size += (points[i+k][l+m].size() - 3)/2;
								count++;
							}
						}
						break;
					}
				}
			}
			
			// find points below latitude
			for( int k = 1; k <= ptRadius && i - k >= 0; k++ ) {
				// find where points are for a given latitude
				for( int l = 1; l < points[i-k].size(); l++ ) {
					if( points[i-k][l][1] > points[i][j][1] || points[i-k].size() - 1 == l ) {
						
						// find points left of point on a given latitude
						for( int m = 0; m <= ptRadius + 1 && l - m >= 0; m++ ) {
							if( Haversine(points[i][j], points[i-k][l-m]) <= epsilon ) {
								distances[i][j].push_back(i-k);
								distances[i][j].push_back(l-m);
								size += (points[i-k][l-m].size() - 3)/2;
								count++;
							}
						}
						
						// find points right of point on a given latitude
						for( int m = 1; m <= ptRadius && points[i-k].size() - l - m > 0; m++ ) {
							if( Haversine(points[i][j], points[i-k][l+m]) <= epsilon ) {
								distances[i][j].push_back(i-k);
								distances[i][j].push_back(l+m);
								size += (points[i-k][l+m].size() - 3)/2;
								count++;
							}
						}
						
						break;
					}
				}
			}
			
			// report size for this point
			distances[i][j].push_back(floor(size/count));
			
		}
	}
	
	
	
	
	
	
	//***** dbscan algorithm *****
	
	// dbscan
	for( unsigned int i = 0; i < points.size(); i++ ) {
		for( unsigned int j = 0; j < points[i].size(); j++ ) {
			// find unlabeled points
			if( points[i][j][2] == 0 ) {
				// if it is a center, create a cluster around it
				if( distances[i][j][distances[i][j].size()-1] >= minPts ) {
					cluster++;
					points[i][j][2] = cluster;
					NbhdJump( i, j, cluster );
				}
				// otherwise mark it as an outlier
				else {
					points[i][j][2] = -1;
				}
			}
		}
	}
	
	
	
	// output results
	cout << cluster << endl;
	
	/*for( unsigned int i = 0; i < points.size(); i++ ) {
		for( unsigned int j = 0; j < points[i].size(); j++ ) {
			cout << points[i][j][0] << "," << points[i][j][1] << "," << points[i][j][2] << "," << distances[i][j][distances[i][j].size() - 1] << endl;
		}
	}*/
	
	out.open("results.geojson");	
	out << "{ \"type\": \"FeatureCollection\",\n\t\"features\": [\n";
	for( int k = 1; k <= cluster; k++ ) {
		out << "\t\t{ \"type\": \"Feature\",\n";
		out << "\t\t\t\"geometry\": {\"type\": \"MultiPoint\", \"coordinates\": [ ";
		bool first = true;
		for( int i = 0; i < points.size(); i++ ) {
			for( unsigned int j = 0; j < points[i].size(); j++ ) {
				if( points[i][j][2] == k ) {
					if( first != true ) {out << ", ";}
					first = false;
					out << "[" << points[i][j][1] << ", " << points[i][j][0] << "]";
				}
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
	
	
	
	
	return 0;
}

void NbhdJump( int i, int j, int cluster ) {
	int size = 0;
	for( int k = 0; k < distances[i][j].size() - 1; k += 2 ) {
		// for each point close to [i][j], find unlabelled points
		if(points[distances[i][j][k]][distances[i][j][k+1]][2] == 0) {
			// for each unlabeled point, if it is dense enough, add it to the cluster and expand
			if(distances[distances[i][j][k]][distances[i][j][k+1]][distances[distances[i][j][k]][distances[i][j][k+1]].size() - 1] >= minPts) {
				points[distances[i][j][k]][distances[i][j][k+1]][2] = cluster;
				NbhdJump(distances[i][j][k], distances[i][j][k+1], cluster);
			}
			// otherwise, mark it as an outlier
			else {
				points[distances[i][j][k]][distances[i][j][k+1]][2] = -1;
			}
		}
	}
}




double Haversine( vector<double> p1, vector<double> p2 ) {
	const double PI = 3.141592653589793;
	double lon1, lon2, lat1, lat2;
	lon1 = (p1[1]*PI)/180;
	lon2 = (p2[1]*PI)/180;
	lat1 = (p1[0]*PI)/180;
	lat2 = (p2[0]*PI)/180;
	
	
	const double R = 6366;
	double a, c, d;
	a = pow( sin((lat1 - lat2)/2), 2 ) + cos(lat1)*cos(lat2)*pow( sin((lon1 - lon2)/2), 2 );
	c = 2 * atan2(sqrt(a), sqrt(1-a));
	d = R * c;
	return d;
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
	
	
	return val;
						
}































































