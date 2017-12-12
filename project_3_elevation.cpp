//
//  project_3_elevation.cpp
//  newtestproj
//  ECE 1305 fall 17
//  Created by David Hoefs on 12/8/17.
//  Copyright © 2017 David Hoefs. All rights reserved.
//  to run this program, set the "string base folder" to the file path which holds the egm files. this program only works with e1 egms, not binary egms.

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cstring>

using namespace std;
//...................Structs...............................
struct Vec3d
{
    double x=0, y=0, z=0;
    
    Vec3d(double my_x=0,double my_y=0,double my_z=0)
    {
        my_x = x;
        my_y=y;
        my_z=z;
    }
};

//...................Function Prototypes...................
void get_filenames(string &fname, string &outname, double &lat, double &lon, double &dlat, double &dlon);
void read_egm(string fname,string &ftype, int &nc, int &nr,float &max, float &min,short * &data_e4, float &Yscale,float &Xscale,float * &data_e1,double* &dnew);
void rescale_EGM_to_PGM(string fname,string outname,string ftype, int &nc, int &nr,float &max, float &min,short * &data, float &Yscale,float &Xscale,float* data_e1,int* &newdat,short* &data_e4);
void illuminate_PGM_from_EGM(string ftype,string outname,int &nc, int &nr,float &max, float &min,short * &data_e4, float &Yscale,float &Xscale,Vec3d v,float* &data_e1,Vec3d sun,int* &newdat,double* &dnew);
Vec3d get_sun_angle();
 Vec3d unit_vector(Vec3d v);
double mag(Vec3d v);
Vec3d cross(Vec3d A, Vec3d B);
double dot_prod(Vec3d A,Vec3d B,double &dot_product);
void write_pgm(string outname,int nc,int nr,float max,float min, float * &data_e1,int* newdat,double* &dnew);
double deg2rad(double deg);
void illum(string outname,int nc,int nr,float max,float min,float Yscale,float Xscale,Vec3d v,double* dnew,int* &newdat);
void write_bin(string outname, int nc, int nr, int* newdat);


 string base_folder = "/Users/DavidHoefs/Desktop/dems_for_project_3/";
int main(void)
{
    string fname,outname,resp;  //declaring my variables.
    string ftype;
    double lat=0,lon=0,dlat=0,dlon=0;
    int nc=0,nr=0;
    float min=0,max=0,Yscale=0,Xscale=0;
    float* data_e1=nullptr;
    short* data_e4=nullptr;
    int* newdat = nullptr;
    double* dnew =nullptr;
    
    Vec3d v,sun;    // declaring some structs
    
    ifstream ifs;       //opening the input filestream
    get_filenames(fname,outname,lat,lon,dlat,dlon);
   
    read_egm(fname,ftype,nc,nr,max,min,data_e4,Yscale,Xscale,data_e1,dnew);
   
    cout << ftype << endl;
    cout << "simple rescale? (y/n): ";      // prompts the user to choose if they want to rescale or illuminate
    cin >> resp;
    if (resp == "y" or resp == "Y")
    rescale_EGM_to_PGM(fname,outname,ftype,nc,nr,max,min,data_e4,Yscale,Xscale,data_e1,newdat,data_e4);
    else if(ftype == "1")
    {
    illuminate_PGM_from_EGM(ftype,outname,nc, nr, max, min, data_e4, Yscale, Xscale, v, data_e1, sun,newdat, dnew);
    write_pgm(outname,nc,nr,max,min,data_e1,newdat,dnew);
    }
    else if(ftype == "4")
    {
        illum(outname,nc,nr,max,min,Yscale,Xscale,v,dnew,newdat);
        write_bin(outname, nc, nr, newdat);
    }
    delete [] data_e1; //delete the memory when done
    delete [] data_e4;
    delete [] newdat;
    delete [] dnew;
   system(("open "+outname+".pgm").c_str());
    return 0;
}
//.......................Functions...............................
void get_filenames(string & fname, string &outname, double &lat, double &lon, double &dlat, double &dlon)
{
    string prompt="enter EGM filename, with NO extension (or q to quit): ";
    string qname;
    bool ok = false;
    
    while (!ok)
    {
        cout << prompt;
        cin >> qname;
        // check to see if we can open the file
        
        fname = base_folder + qname + ".egm";
        ifstream ifs(fname);
        if (!ifs)
        {
            if (qname == "q" || qname == "Q")
                exit(1);        // bail out if user wants to quit.
            prompt = "error opening " + fname + "\nenter EGM filename, with NO extension (or q to quit): ";
        }
        else
        {
            ok = true;      // ok. we're good.  we can open the file for reading.
            ifs.close();
        }
    }
    
       cout << "enter output PGM filename: ";
        cin >> outname;
    outname = base_folder + outname;
    //outname = base_folder + qname + "_illum.pgm";
    
    if (!strncmp(qname.c_str(), "mola",4))
    {
        cout << "lat, lon (deg)? ";
        cin >> lat >> lon;
        cout << "dlat, dlon (deg)? ";
        cin >> dlat >> dlon;
    }
}
//............................Function to read the egm.................................
// this function reads the egm and determines the header information, then fills an array with the data in the file.
void read_egm(string fname,string &ftype, int &nc, int &nr,float &max, float &min,short * &data_e4, float&Yscale,float &Xscale,float * &data_e1,double * &dnew)
{
    ifstream ifs;           //opening the input filestream
    ifs.open(fname);        // open the file
    if(!ifs)                // checks to make sure the file was successfully opened
    {
        cerr << "cannot open file" << endl;
    }
    
    char* type = new char[2];       // declaring a char pointer to determine if the file is e1 or e4
    ifs.read(type,2);
    ftype = type[1];                // assigns the variable ftype to the value of type (e1 or e4)
    cout << ftype << endl;
    delete [] type;                 // cleaning up when Im finished with the char pointer.
    type = nullptr;                 // setting type to null
    
        if(ftype == "1")        // if the file was e1
        {
        ifs >> nr >> nc >> min >> max >> Xscale >>Yscale;   // read in the header information using >> operator
        float size = nc*nr;                         // declaring the variable size to be equal to total number of values in file
        data_e1 = new float[size];      // creating a pointer array to hold the data in the file.
        cout << nc << " " << nr << endl;;
        for(int i=0;i<size;i++)         // for loop to fill the dynamically allocated array with values.
        {
            
            ifs >> data_e1[i];
            
        }
        cout << data_e1[0] << endl;
        }else if(ftype == "4")
        {
            ifs.close();
            ifs.open(fname,ios::binary);
    type = new char[8];
    ifs.read(type,8);
    cout << type << endl;
    int * header = new int[6];
    ifs.read(reinterpret_cast<char*>(header),6*sizeof(int));
    
    nc = header[0];
    nr = header[1];
    min = header[2];
    max = header[3];
    Xscale= header[4];
    Yscale = header[5];
    cout << nc << " " << nr << " " << min << " " <<max << " " << Xscale << " " <<Yscale << endl;
    int size = nc*nr;
    data_e4 = new short[size];
    ifs.read(reinterpret_cast<char*>(data_e4),size * sizeof(short));
    dnew = new double[size];
    for(int i=0;i<size;i++)
    {
        dnew[i] = data_e4[i];
    }
   
    ifs.close();
        }
}
        
    
         // close the input filestream when finished.
    
    

    
        
        
    

//..................................Function to do a simple rescale of the elevation values.........................

void rescale_EGM_to_PGM(string fname,string outname ,string ftype, int &nc, int &nr,float &max, float &min,short * &data, float &Yscale,float &Xscale,float* data_e1,int* &newdat,short* &data_e4)
{
    
    cout << ftype << endl;
    ofstream opfs;                  // declating the output filestream
    if(ftype == "1")                // if the file type was e1
    {
        
    opfs.open(outname+".pgm");        // open file to write to
    opfs << setprecision(0) << fixed; // setting the precision so float values will be whole numbers.
    opfs << "P2" << endl;          // writing out the header information to the file
    opfs << nr << endl;
    opfs << nc << endl;
    opfs << 255 << endl;
    for(int i=0;i<nr;i++)         //for number of rows
    {
        for(int j=0;j<nc;j++)       // for number of columns
        {
            float val = data_e1[i*nc+j];     // finding the value of each pix
            float bness = (val - min)/(max-min); // value - min elev / max elev - min elev multiplied by 255
            opfs << bness*255 << "\t";      // outputting the rescaled brightness to the file
        }opfs << endl;                      // starting new line after each row
    }
    
    
    
    }
}
//...............function to convery degrees to radians...................
// input: degrees
//output: radians
double deg2rad(double deg)

{
    
    return M_PI/180.0 * deg;
    
}

//.............Function to get the angle of the sun in x,y,z coordinates....................
// prompts user for azimuth and elevation of the sun, then converts the user input to (x,y,z) coordinates
// input: azimuth and elevation
// output: returns a Vec3d with x,y,z coordinates
Vec3d get_sun_angle(void)
{
    // get direction to sun
    double az, el;
    cout << "sun az (deg): ";
    cin >> az;
     az = deg2rad(az);
    cout << "sun el (deg): ";
    cin >> el;
    el = deg2rad(el);
    
    // convert it to a 3-D vector:
    Vec3d sun;
    sun.x = cos(el) * cos(az);
    sun.y = cos(el) * sin(az);
    sun.z = sin(el);
    return sun;
}
//.......................Function to find unit vector (Normalized vector)......................
// input: a Vec3d
// output: the normalized version of the Vec3d you entered
Vec3d unit_vector(Vec3d v)
{
    double mymag = mag(v); // mymag set equal to the magnitude
    v.x /= mymag;           // divides each coordinate by its magnitude
    v.y /= mymag;
    v.z /= mymag;
    return v;
}
//...................Function to find the magnitude of a vector........................
//input: a Vec3d vector
//output: a double that is the magnitude
double mag(Vec3d v)
{
    return sqrt(v.x*v.x + v.y*v.y + v.z * v.z); // magnitude is the square root of each coordinate squared
}
//.......................Function to illuminated egm values........................
// input: nc,nr,max,min, short pointer array, yscale,xscale, Vec3d v and sun, and float pointer array;
// this function finds two vectors adjacent to each other, calculates their surface normal, then finds the cos of the angle between the surface normal and the ray coming from the sun
void illuminate_PGM_from_EGM(string ftype,string outname,int &nc, int &nr,float &max, float &min,short * &data_e4, float &Yscale,float &Xscale,Vec3d v,float* &data_e1,Vec3d sun, int* &newdat,double* &dnew)
{
    Vec3d a,b;
    sun = get_sun_angle();   // prompts user for sun azimuth and elevation
    
    int* prev = new int[nc*nr];
    
    newdat = new int[nc*nr];
    for(int r =0;r<nr;r++)      // for each row
    {
        for(int c=0;c<nc;c++)// for each column
        {
            a.x = Xscale;       // setting x in vec a to be the xscaling
            a.y=1;              // setting y in vec a to be 1;
            int pointa = data_e1[r*nc+c+1] - data_e1[r*nc+c];    // finding the difference of a value and the value to its right
            a.z= pointa;        // setting z in vec a to the difference of the points
            b.x=1;              // setting x in vec b to 1
            b.y=Yscale;         // setting y in vec b to the yscaling
            int pointb = data_e1[r+1*nc+c] - data_e1[r*nc+c]; // finding the change in elevation from the original elevation and the elevation directly beneath it
            b.z = pointb;  // setting z in vec b to that value that was found above
            
            Vec3d anorm = unit_vector(a);   // finding the unit vector of vec a
            Vec3d bnorm = unit_vector(b);   // finding the unit vector of vec b
            Vec3d snorm = cross(bnorm,anorm); // finding the cross product of normalized vectors b and a
            double brightness = dot_prod(snorm,sun,brightness); // finding the brightness at each point by calculating the dot product of the surface normal and the normalized vector of the sun
            
            float bness = ((.9*brightness)+.1) *255;    // bness(brightness value) is the dot product times 255 (to get a value between 0 and 255 for pgm) .9 and +.1 added to represent the sky brighntess
         
           // if the brighness is negative (cosine of the angle was negative, meaning the brighness of that point is in a shadow) set the brightness to a positive number
                if(bness <0)
                {
                    
                    bness = bness + 100 ;
                    
                }
            
            
            newdat[r*nc + c] = bness;
            prev[r*nc +c] = bness;
            
        }
      
    }
    
   
    delete [] prev;
    prev = nullptr;
    }
//....................... Function to find the cross product of two vectors........................
//input: two Vec3d vectors
// output: a new Vec3d vector
Vec3d cross(Vec3d A,Vec3d B)
{
    Vec3d snorm;
    snorm.x=(A.y*B.z)-(B.y*A.z);
    snorm.y=(A.z*B.x) - (B.z*A.x);
    snorm.z=(A.x*B.y) - (B.x*A.y);
    return snorm;
}
//...................Function to find dot product of two vectors....................
//input: two vectors
//output: a scalar value that is the dot product
double dot_prod(Vec3d A,Vec3d B,double &dot_product)
{
    dot_product = (A.x*B.x) + (A.y*B.y)+ (A.z*B.z);
    return dot_product;
    
}
void write_pgm(string outname,int nc,int nr,float max,float min, float * &data_e1,int* newdat, double* &dnew)
{
    
    int j=0;
    ofstream opfs;
    opfs.open(outname+".pgm");
    opfs << "P2" << endl;
    opfs << nr << endl;
    opfs << nc << endl;
    opfs << 255 << endl;
  
    for(int r=0;r<nr;r++)
    {
        for(int c=0;c<nc;c++)
        {
         
         
            opfs << newdat[r*nc + c];
            
            opfs << "\t";
        }opfs << endl;
    }
    opfs.close();
}
void illum(string outname,int nc,int nr,float max,float min,float Yscale,float Xscale, Vec3d v,double* dnew,int* &newdat)
{
    Vec3d a,b,sun;
    sun = get_sun_angle();   // prompts user for sun azimuth and elevation
    
    newdat = new int[nc*nr];
    
  
    for(int r =0;r<nr-1;r++)      // for each row
    {
        for(int c=0;c<nc-1;c++)// for each column
        {
            a.x = Xscale;       // setting x in vec a to be the xscaling
            a.y=1;              // setting y in vec a to be 1;
           double pointa = dnew[r*nc+c+1] - dnew[r*nc+c];    // finding the difference of a value and the value to its right
            a.z= pointa;        // setting z in vec a to the difference of the points
            b.x=1;              // setting x in vec b to 1
            b.y=Yscale;         // setting y in vec b to the yscaling
            double pointb = dnew[r+1*nc+c] - dnew[r*nc+c]; // finding the change in elevation from the original elevation and the elevation directly beneath it
            b.z = pointb;  // setting z in vec b to that value that was found above
            
            Vec3d anorm = unit_vector(a);   // finding the unit vector of vec a
            Vec3d bnorm = unit_vector(b);   // finding the unit vector of vec b
            Vec3d snorm = cross(bnorm,anorm); // finding the cross product of normalized vectors b and a
            double brightness = dot_prod(snorm,sun,brightness); // finding the brightness at each point by calculating the dot product of the surface normal and the normalized vector of the sun
            
            int bness = ((.9*brightness)+.1) *255;    // bness(brightness value) is the dot product times 255 (to get a value between 0 and 255 for pgm) .9 and +.1 added to represent the sky brighntess
            
            // if the brighness is negative (cosine of the angle was negative, meaning the brighness of that point is in a shadow) set the brightness to a positive number
            if(bness <0)
            {
                
                bness = 115;
                
            }
            
            
            newdat[r*nc + c] = bness;
          
            
        }
        
    }
    
  
}
void write_bin(string outname, int nc, int nr, int* newdat)
{
    ofstream opfs;
    unsigned char* cdata = new unsigned char[nc*nr];
    for(int i=0;i<nc*nr;i++)
    {
        cdata[i]=newdat[i];
    }
    opfs.open(outname+".pgm",ios::binary);
    opfs << "P5" << endl;
    opfs << nr << endl;
    opfs << nc << endl;
    opfs << 255 << endl;
    opfs.write(reinterpret_cast<char*>(cdata),nc*nr*sizeof(unsigned char));
    opfs.close();
}
    

