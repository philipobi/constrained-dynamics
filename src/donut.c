#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#define x_stretch 1.5
#define screeny 50

double R=2.0;
double r=1.0;
double d1=20.0;
double d2=4.0;
int n_theta = 100;
int n_phi = 100;

int i,j,y_index,x_index;

const int screenx = screeny*x_stretch;
int half_screen_x = screenx/2;
int half_screen_y = screeny/2;

double cos_psi,one_min_cos_psi,sin_psi,sin_phi,sin_theta,cos_phi,cos_theta,x,y,z,factor,theta,phi,psi,dtheta,dphi,dpsi,n1,n2,n3,Rrcosphi;

char screen[screeny][screenx+1];
char screen_empty[screeny][screenx+1];
double xbuf[screeny][screenx];
double xbuf_empty[screeny][screenx];

int light_index;

double rot_matrix[3][3];
double r_vec[3];
double r_vec_rot[3];
double normal_vec[3];
double normal_vec_rot[3];
double light_direction[3];

const char light_level[] = {'.', ',', '-', '~', ':', ';', '=', '!', '*', '#', '$', '@'};


void setup()
{
    for(i=0;i<screeny;++i)
    for(j=0;j<screenx;++j)
    {
        screen_empty[i][j]=' ';
        if(j==screenx-1) screen_empty[i][screenx]='\0';
    }
}

double dot(const double v1[3], const double v2[3])
{
    double sum=0.0;
    for(int k=0;k<3;++k) sum += v1[k]*v2[k];
    return sum;
}

void mat_vec_mul(const double mat[3][3], const double vec[3], double* result)
{
    double sum;
    for(int i=0; i<3; ++i)
    {
        sum = 0.0;
        for(int j=0; j<3; ++j) sum+=mat[i][j]*vec[j];
        result[i]=sum;
    }
}

void render(double psi)
{
    memcpy(screen,screen_empty,sizeof(char)*screeny*(screenx+1));
    memcpy(xbuf,xbuf_empty,sizeof(double)*screeny*screenx);
    
    sin_psi = sin(psi);
    cos_psi = cos(psi);
    one_min_cos_psi = 1.0-cos_psi;

    rot_matrix[0][0] = cos_psi+n1*n1*one_min_cos_psi;
    rot_matrix[0][1] = n1*n2*one_min_cos_psi-n3*sin_psi;
    rot_matrix[0][2] = n1*n3*one_min_cos_psi+n2*sin_psi;
    
    rot_matrix[1][0] = n2*n1*one_min_cos_psi+n3*sin_psi;
    rot_matrix[1][1] = cos_psi+n2*n2*one_min_cos_psi;;
    rot_matrix[1][2] = n2*n3*one_min_cos_psi-n1*sin_psi;

    rot_matrix[2][0] = n3*n1*one_min_cos_psi-n2*sin_psi;
    rot_matrix[2][1] = n3*n2*one_min_cos_psi+n1*sin_psi;
    rot_matrix[2][2] = cos_psi+n3*n3*one_min_cos_psi;

    theta=0.0; phi=0.0;

    for(i=0;i<n_theta;++i,theta+=dtheta)
    {
        sin_theta = sin(theta);
        cos_theta = cos(theta);

        for(j=0;j<n_phi;++j,phi+=dphi)
        {
            sin_phi=sin(phi);
            cos_phi=cos(phi);
            
            Rrcosphi = R+r*cos_phi;
            
            r_vec[0]=-r*sin_phi;
            r_vec[1]=Rrcosphi*cos_theta;
            r_vec[2]=Rrcosphi*sin_theta;

            mat_vec_mul(rot_matrix,r_vec,r_vec_rot);
            
            x = -r_vec_rot[0]+d2;
            y = r_vec_rot[1];
            z = r_vec_rot[2];
            factor = d1/x;
            
            x_index = half_screen_x + round(y*x_stretch*factor);
            y_index = half_screen_y + round(z*factor);

            if(factor > xbuf[y_index][x_index])
            {
                normal_vec[0] = -sin_phi;
                normal_vec[1] = cos_phi*cos_theta;
                normal_vec[2] = cos_phi*sin_theta;
                mat_vec_mul(rot_matrix,normal_vec,normal_vec_rot);
                
                xbuf[y_index][x_index] = factor;
                screen[y_index][x_index] = light_level[(int) round((dot(normal_vec_rot,light_direction)*(-1.0)+1.0)*5.5)];
            }
        }
    }
    for(int i=0; i<screeny; ++i) printf("%s\n",screen[i]);
}

int main()
{   
    setup();
    
    n1=0.0;n2=0.0;n3=1.0;
    dtheta = 2*M_PI/n_theta;
    dphi = 2*M_PI/n_phi;
    dpsi = 0.001*M_PI;

    light_direction[0]=0.;
    light_direction[1]=-1.;
    light_direction[2]=0.;

    for(psi=0.0;;psi+=dpsi)
    {
        render(psi);
        usleep(1000);
        printf("\033[%dA",screeny);
    }
    
}