#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define grid 51

int main()
{
    double u[grid][grid+1],un[grid][grid+1],du[grid][grid+1],u0[grid][grid+1],u00[grid][grid+1];
    double v[grid+1][grid],vn[grid+1][grid],dv[grid+1][grid],v0[grid+1][grid],v00[grid+1][grid];
    double p[grid+1][grid+1],pc[grid+1][grid+1],pc0[grid+1][grid+1],pcn[grid+1][grid+1],bc[grid+1][grid+1],A[grid+1][grid+1];
    double ap[grid+1][grid+1],ae[grid+1][grid+1],aw[grid+1][grid+1],an[grid+1][grid+1],as[grid+1][grid+1],B[grid+1][grid+1];
    double m[grid+1][grid+1],a[grid+1],b[grid+1],c[grid+1],d[grid+1],gamma[grid+1],beta[grid+1];
    double Re = 1000.0;
    double error = 1.0, errorT = 1.0;
    double dx = 1.0/(grid-1), dy=1.0/(grid-1), dt = 0.01;
    int i,j,step = 1;

    // initialization of u
    for(i=1;i<=grid-2;i++)
    {
        for(j=0;j<=grid-2;j++)
            u0[i][j]=0.0;
        u0[i][grid]=1.0;
        u0[i][grid-1]=1.0;
    }
    for(j=0;j<=grid;j++)
    {
        u0[0][j] = 0.0;
        u0[grid-1][j] = 0.0;
    }

    // initialization of v
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid-1;j++)
            v0[i][j]=0.0;
    }

    //initialization of p
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid;j++)
            p[i][j]=1.0;
    }

    //initialization of pc
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid;j++)
            pc[i][j]=0.0;
    }
    //while(step<=100)
    //{
        for(i=0;i<=grid-1;i++)
        {
            for(j=0;j<=grid;j++)
                u[i][j]=u0[i][j];
        }

        for(i=0;i<=grid;i++)
        {
            for(j=0;j<=grid-1;j++)
                v[i][j]=v0[i][j];
        }
        error = 1.0;
        while(error>0.0000001)
        {
            for(i=0;i<=grid-1;i++)
            {
                for(j=0;j<=grid;j++)
                    u[i][j]=u0[i][j];
            }

            for(i=0;i<=grid;i++)
            {
                for(j=0;j<=grid-1;j++)
                    v[i][j]=v0[i][j];
            }
            //x-momentum equation
            errorT = 1.0;
            while(errorT>0.0000001)
            {
                //x-momentum equation coefficients
                for(i=1;i<=grid-2;i++)
                {
                    for(j=1;j<=grid-1;j++)
                    {
                        ap[i][j] =  (1.0*dx*dy/dt)
                                    +(u[i+1][j]+u[i][j])*dy/(4.0) - (u[i-1][j]+u[i][j])*dy/(4.0)
                                    +(v[i+1][j]+v[i][j])*dx/(4.0) - (v[i+1][j-1]+v[i][j-1])*dx/(4.0)
                                    + 2.0*dy/(Re*dx) + 2.0*dx/(Re*dy);
                        ae[i][j] = -(u[i+1][j]+u[i][j])*dy/(4.0) + 1.0*dy/(Re*dx);
                        aw[i][j] = (u[i-1][j]+u[i][j])*dy/(4.0) + 1.0*dy/(Re*dx);
                        an[i][j] = -(v[i+1][j]+v[i][j])*dx/(4.0) + 1.0*dx/(Re*dy);
                        as[i][j] = (v[i+1][j-1]+v[i][j-1])*dx/(4.0) + 1.0*dx/(Re*dy);
                        A[i][j] = 1.0*dy;
                        du[i][j] = A[i][j]/ap[i][j];
                        B[i][j] = u[i][j]*dx*dy/dt;
                    }
                }
                //solving x-momentum equation
                for(j=1;j<=grid-1;j++)
                {
                    b[1] = ap[1][j];
                    c[1] = -ae[1][j];
                    d[1] = an[1][j]*u0[1][j+1] + as[1][j]*u0[1][j-1] + A[1][j]*(p[1][j] - p[2][j]) + B[i][j];
                    for(i=2;i<=grid-3;i++)
                    {
                        a[i] = -aw[i][j];
                        b[i] = ap[i][j];
                        c[i] = -ae[i][j];
                        d[i] = an[i][j]*u0[i][j+1] + as[i][j]*u0[i][j-1] + A[i][j]*(p[i][j] - p[i+1][j]) + B[i][j];
                    }
                    a[grid-2] = -aw[grid-2][j];
                    b[grid-2] = ap[grid-2][j];
                    d[grid-2] = an[grid-2][j]*u0[grid-2][j+1] + as[grid-2][j]*u0[grid-2][j-1] + A[grid-2][j]*(p[grid-2][j] - p[grid-1][j]) + B[i][j];
                    gamma[1] = c[1]/b[1];
                    beta[1] = d[1]/b[1];
                    for(i=2;i<=grid-2;i++)
                    {
                        if(i!=grid-2)
                            gamma[i] = c[i]/(b[i]-a[i]*gamma[i-1]);
                        beta[i] = (d[i]-a[i]*beta[i-1])/(b[i]-a[i]*gamma[i-1]);
                    }
                    un[grid-2][j] = beta[grid-2];

                    for(i=grid-3;i>=1;i--)
                        un[i][j] = beta[i] - gamma[i]*un[i+1][j];
                }
                for(j=0;j<=grid;j++)
                {
                    un[0][j] = 0.0;
                    un[grid-1][j] = 0.0;
                }
                for(i=1;i<=grid-2;i++)
                {
                    un[i][0] = -un[i][1];
                    un[i][grid] = 2.0-un[i][grid-1];
                }

                //error calculated for tdma iteration
                errorT=0.0;
                for(i=1;i<=grid-2;i++)
                {
                    for(j=1;j<=grid-1;j++)
                    {
                        m[i][j] = ( un[i][j]-u0[i][j] );
                        errorT = errorT + fabs(m[i][j]);
                    }
                }
                //next tdma iteration
                for(i=0;i<=grid-1;i++)
                {
                    for(j=0;j<=grid;j++)
                        u0[i][j]=(u0[i][j]+un[i][j])/2.0;
                }

            }

            //y-momentum equation
            errorT = 1.0;
            while(errorT>0.0000001)
            {
                //y-momentum equation coefficients
                for(i=1;i<=grid-1;i++)
                {
                    for(j=1;j<=grid-2;j++)
                    {
                        ap[i][j] =  (1.0*dx*dy/dt)
                                    +(u[i][j+1]+u[i][j])*dy/(4.0) - (u[i-1][j+1]+u[i-1][j])*dy/(4.0)
                                    +(v[i][j+1]+v[i][j])*dy/(4.0) - (v[i][j-1]+v[i][j])*dy/(4.0)
                                    + 2.0*dy/(Re*dx) + 2.0*dx/(Re*dy);
                        ae[i][j] = -(u[i][j+1]+u[i][j])*dy/(4.0) + 1.0*dy/(Re*dx);
                        aw[i][j] = (u[i-1][j+1]+u[i-1][j])*dy/(4.0) + 1.0*dy/(Re*dx);
                        an[i][j] = -(v[i][j+1]+v[i][j])*dx/(4.0) + 1.0*dx/(Re*dy);
                        as[i][j] = (v[i][j-1]+v[i][j])*dx/(4.0) + 1.0*dx/(Re*dy);
                        A[i][j] = 1.0*dx;
                        dv[i][j] = A[i][j]/ap[i][j];
                        B[i][j] = v[i][j]*dx*dy/dt;
                    }
                }

                //solving y-momentum equation
                for(i=1;i<=grid-1;i++)
                {
                    b[1] = ap[i][1];
                    c[1] = -an[i][1];
                    d[1] = ae[i][1]*v0[i+1][1] + aw[i][1]*v0[i-1][1] + A[i][1]*(p[i][1] - p[i][2]) + B[i][j];
                    for(j=2;j<=grid-3;j++)
                    {
                        a[j] = -as[i][j];
                        b[j] = ap[i][j];
                        c[j] = -an[i][j];
                        d[j] = ae[i][j]*v0[i+1][j] + aw[i][j]*v0[i-1][j] + A[i][j]*(p[i][j] - p[i][j+1]) + B[i][j];
                    }
                    a[grid-2] = -as[i][grid-2];
                    b[grid-2] = ap[i][grid-2];
                    d[grid-2] = ae[i][grid-2]*v0[i+1][grid-2] + aw[i][grid-2]*v0[i-1][grid-2] + A[i][grid-2]*(p[i][grid-2] - p[i][grid-1]) + B[i][j];
                    gamma[1] = c[1]/b[1];
                    beta[1] = d[1]/b[1];
                    for(j=2;j<=grid-2;j++)
                    {
                        if(j!=grid-2)
                            gamma[j] = c[j]/(b[j]-a[j]*gamma[j-1]);
                        beta[j] = (d[j]-a[j]*beta[j-1])/(b[j]-a[j]*gamma[j-1]);
                    }
                    vn[i][grid-2] = beta[grid-2];
                    for(j=grid-3;j>=1;j--)
                        vn[i][j] = beta[j] - gamma[j]*vn[i][j+1];
                }
                for(i=0;i<=grid;i++)
                {
                    vn[i][0] = 0.0;
                    vn[i][grid-1] = 0.0;
                }
                for(j=1;j<=grid-2;j++)
                {
                    vn[0][j] = -vn[1][j];
                    vn[grid][j] = -vn[grid-1][j];
                }

                //error calculated for tdma iteration
                errorT=0.0;
                for(i=1;i<=grid-1;i++)
                {
                    for(j=1;j<=grid-2;j++)
                    {
                        m[i][j] = ( vn[i][j]-v0[i][j] );
                        errorT = errorT + fabs(m[i][j]);
                    }
                }
                //next tdma iteration
                for(i=0;i<=grid;i++)
                {
                    for(j=0;j<=grid-1;j++)
                        v0[i][j]=(v0[i][j]+vn[i][j])/2.0;
                }

            }

            //pressure correction equation

            //initialization of pressure correction
            for(i=1;i<=grid-1;i++)
            {
                for(j=1;j<=grid-1;j++)
                    pc[i][j]=0.0;
            }

            errorT = 1.0;
            while(errorT>0.0000001)
            {

                //pressure correction coefficients

                //coefficients for all nodes except for left, right, top and bottom wall
                for(i=2;i<=grid-2;i++)
                {
                    for(j=2;j<=grid-2;j++)
                    {
                        ap[i][j] = du[i][j]+du[i-1][j]+dv[i][j]+dv[i][j-1];
                        ae[i][j] = du[i][j];
                        aw[i][j] = du[i-1][j];
                        an[i][j] = dv[i][j];
                        as[i][j] = dv[i][j-1];
                        bc[i][j] = un[i-1][j]-un[i][j]+vn[i][j-1]-vn[i][j];
                    }
                }
                //coefficients for left and right wall
                for(j=2;j<=grid-2;j++)
                {
                    ap[1][j] = du[1][j]+dv[1][j]+dv[1][j-1];
                    ae[1][j] = du[1][j];
                    aw[1][j] = 0.0;
                    an[1][j] = dv[1][j];
                    as[1][j] = dv[1][j-1];
                    bc[1][j] = u0[0][j]-un[1][j]+vn[1][j-1]-vn[1][j];
                    ap[grid-1][j] = du[grid-2][j]+dv[grid-1][j]+dv[grid-1][j-1];
                    ae[grid-1][j] = 0.0;
                    aw[grid-1][j] = du[grid-2][j];
                    an[grid-1][j] = dv[grid-1][j];
                    as[grid-1][j] = dv[grid-1][j-1];
                    bc[grid-1][j] = un[grid-2][j]-u0[grid-1][j]+vn[grid-1][j-1]-vn[grid-1][j];
                }
                //coefficients for top and bottom wall
                for(i=2;i<=grid-2;i++)
                {
                    ap[i][1] = du[i][1]+du[i-1][1]+dv[i][1];
                    ae[i][1] = du[i][1];
                    aw[i][1] = du[i-1][1];
                    an[i][1] = dv[i][1];
                    as[i][1] = 0.0;
                    bc[i][1] = un[i-1][1]-un[i][1]+v0[i][0]-vn[i][1];
                    ap[i][grid-1] = du[i][grid-1]+du[i-1][grid-1]+dv[i][grid-2];
                    ae[i][grid-1] = du[i][grid-1];
                    aw[i][grid-1] = du[i-1][grid-1];
                    an[i][grid-1] = 0.0;
                    as[i][grid-1] = dv[i][grid-2];
                    bc[i][grid-1] = un[i-1][grid-1]-un[i][grid-1]+vn[i][grid-2]-v0[i][grid-1];
                }
                //coefficients for bottom-left node
                ap[1][1] = du[1][1]+dv[1][1];
                ae[1][1] = du[1][1];
                aw[1][1] = 0.0;
                an[1][1] = dv[1][1];
                as[1][1] = 0.0;
                bc[1][1] = u0[0][1]-un[1][1]+v0[1][0]-vn[1][1];
                //coefficients for top-left node
                ap[1][grid-1] = du[1][grid-1]+dv[1][grid-2];
                ae[1][grid-1] = du[1][grid-1];
                aw[1][grid-1] = 0.0;
                an[1][grid-1] = 0.0;
                as[1][grid-1] = dv[1][grid-2];
                bc[1][grid-1] = u0[0][grid-1]-un[1][grid-1]+vn[1][grid-2]-v0[1][grid-1];
                //coefficients for bottom-right node
                ap[grid-1][1] = du[grid-2][1]+dv[grid-1][1];
                ae[grid-1][1] = 0.0;
                aw[grid-1][1] = du[grid-2][1];
                an[grid-1][1] = dv[grid-1][1];
                as[grid-1][1] = 0.0;
                bc[grid-1][1] = un[grid-2][1]-u0[grid-1][1]+v0[grid-1][0]-vn[grid-1][1];
                //coefficients for top-right node
                ap[grid-1][grid-1] = du[grid-2][grid-1]+dv[grid-1][grid-2];
                ae[grid-1][grid-1] = 0.0;
                aw[grid-1][grid-1] = du[grid-2][grid-1];
                an[grid-1][grid-1] = 0.0;
                as[grid-1][grid-1] = dv[grid-1][grid-2];
                bc[grid-1][grid-1] = un[grid-2][grid-1]-u0[grid-1][grid-1]+vn[grid-1][grid-2]-v0[grid-1][grid-1];

                //solving pressure equation
                for(i=1;i<=grid-1;i++)
                {
                    //solving for left wall
                    if(i==1)
                    {
                    b[1] = ap[i][1];
                    c[1] = -an[i][1];
                    d[1] = ae[i][1]*pc[i+1][1] + bc[i][1];
                    for(j=2;j<=grid-2;j++)
                    {
                        a[j] = -as[i][j];
                        b[j] = ap[i][j];
                        c[j] = -an[i][j];
                        d[j] = ae[i][j]*pc[i+1][j] + bc[i][j];
                    }
                    a[grid-1] = -as[i][grid-1];
                    b[grid-1] = ap[i][grid-1];
                    d[grid-1] = ae[i][grid-1]*pc[i+1][grid-1] + bc[i][grid-1];
                    gamma[1] = c[1]/b[1];
                    beta[1] = d[1]/b[1];
                    for(j=2;j<=grid-1;j++)
                    {
                        if(j!=grid-1)
                            gamma[j] = c[j]/(b[j]-a[j]*gamma[j-1]);
                        beta[j] = (d[j]-a[j]*beta[j-1])/(b[j]-a[j]*gamma[j-1]);
                    }
                    pcn[i][grid-1] = beta[grid-1];
                    for(j=grid-2;j>=1;j--)
                        pcn[i][j] = beta[j] - gamma[j]*pcn[i][j+1];
                    }

                    //solving for right wall
                    else if(i==grid-1)
                    {
                        b[1] = ap[i][1];
                        c[1] = -an[i][1];
                        d[1] = aw[i][1]*pc[i-1][1] + bc[i][1];
                        for(j=2;j<=grid-2;j++)
                        {
                            a[j] = -as[i][j];
                            b[j] = ap[i][j];
                            c[j] = -an[i][j];
                            d[j] = aw[i][j]*pc[i-1][j] + bc[i][j];
                        }
                        a[grid-1] = -as[i][grid-1];
                        b[grid-1] = ap[i][grid-1];
                        d[grid-1] = aw[i][grid-1]*pc[i-1][grid-1] + bc[i][grid-1];
                        gamma[1] = c[1]/b[1];
                        beta[1] = d[1]/b[1];
                        for(j=2;j<=grid-1;j++)
                        {
                            if(j!=grid-1)
                                gamma[j] = c[j]/(b[j]-a[j]*gamma[j-1]);
                            beta[j] = (d[j]-a[j]*beta[j-1])/(b[j]-a[j]*gamma[j-1]);
                        }
                        pcn[i][grid-1] = beta[grid-1];
                        for(j=grid-2;j>=1;j--)
                            pcn[i][j] = beta[j] - gamma[j]*pcn[i][j+1];
                    }

                    //solving for all nodes except for left and right wall
                    else
                    {
                        b[1] = ap[i][1];
                        c[1] = -an[i][1];
                        d[1] = ae[i][1]*pc[i+1][1] + aw[i][1]*pc[i-1][1] + bc[i][1];
                        for(j=2;j<=grid-2;j++)
                        {
                            a[j] = -as[i][j];
                            b[j] = ap[i][j];
                            c[j] = -an[i][j];
                            d[j] = ae[i][j]*pc[i+1][j] + aw[i][j]*pc[i-1][j] + bc[i][j];
                        }
                        a[grid-1] = -as[i][grid-1];
                        b[grid-1] = ap[i][grid-1];
                        d[grid-1] = ae[i][grid-1]*pc[i+1][grid-1] + aw[i][grid-1]*pc[i-1][grid-1] + bc[i][grid-1];
                        gamma[1] = c[1]/b[1];
                        beta[1] = d[1]/b[1];
                        for(j=2;j<=grid-1;j++)
                        {
                            if(j!=grid-1)
                                gamma[j] = c[j]/(b[j]-a[j]*gamma[j-1]);
                            beta[j] = (d[j]-a[j]*beta[j-1])/(b[j]-a[j]*gamma[j-1]);
                        }
                        pcn[i][grid-1] = beta[grid-1];
                        for(j=grid-2;j>=1;j--)
                            pcn[i][j] = beta[j] - gamma[j]*pcn[i][j+1];
                    }

                }

                //error calculated for tdma iteration
                errorT=0.0;
                for(i=1;i<=grid-1;i++)
                {
                    for(j=1;j<=grid-1;j++)
                    {
                        m[i][j] = ( pcn[i][j]-pc[i][j] );
                        errorT = errorT + fabs(m[i][j]);
                    }
                }
                //next iteration of tdma
                for(i=1;i<=grid-1;i++)
                {
                    for(j=1;j<=grid-1;j++)
                        pc[i][j] = (pc[i][j]+pcn[i][j])/2.0;
                }

            }


            //adding correction to u
            for(i=1;i<=grid-2;i++)
            {
                for(j=1;j<=grid-1;j++)
                    u0[i][j] = un[i][j] + du[i][j]*(pcn[i][j] - pcn[i+1][j]);
            }
            for(j=0;j<=grid;j++)
            {
                u0[0][j] = 0.0;
                u0[grid-1][j] = 0.0;
            }
            for(i=1;i<=grid-2;i++)
            {
                u0[i][0] = -u0[i][1];
                u0[i][grid] = 2.0-u0[i][grid-1];
            }

            //adding correction to v
            for(i=1;i<=grid-1;i++)
            {
                for(j=1;j<=grid-2;j++)
                    v0[i][j] = vn[i][j] + dv[i][j]*(pcn[i][j] - pcn[i][j+1]);
            }
            for(i=0;i<=grid;i++)
            {
                v0[i][0] = 0.0;
                v0[i][grid-1] = 0.0;
            }
            for(j=1;j<=grid-2;j++)
            {
                v0[0][j] = -v0[1][j];
                v0[grid][j] = -v0[grid-1][j];
            }

            //error
            error = 0.0;
            for(i=1;i<=grid-2;i++)
            {
                for(j=1;j<=grid-2;j++)
                {
                    m[i][j] = ( (un[i-1][j]-un[i][j])*dy
                               +(vn[i][j-1]-vn[i][j])*dx );
                    error = error + fabs(m[i][j]);
                }
            }
            printf("Error is %5.9lf for step %d\n",error,step);

            //adding correction to p
            for(i=1;i<=grid-1;i++)
            {
                for(j=1;j<=grid-1;j++)
                    p[i][j] = p[i][j] + 0.01*pcn[i][j];
            }
            for(j=1;j<=grid-1;j++)
            {
                p[0][j] = p[1][j];
                p[grid][j] = p[grid-1][j];
            }
            for(i=0;i<=grid;i++)
            {
                p[i][0] = p[i][1];
                p[i][grid] = p[i][grid-1];
            }
            step++;
        }
        double u_c[grid][grid], v_c[grid][grid], p_c[grid][grid];

        for(i=0;i<=grid-1;i++)
        {
            for(j=0;j<=grid-1;j++)
            {
                u_c[i][j] = 0.5*(u0[i][j]+u0[i][j+1]);
                v_c[i][j] = 0.5*(v0[i][j]+v0[i+1][j]);
                p_c[i][j] = 0.25*(p[i][j]+p[i][j+1]+p[i+1][j]+p[i+1][j+1]);
            }
        }


        // OUTPUT DATA
        FILE *fout2, *fout3;
        fout2 = fopen("UVP_steady.plt","a+t");
        fout3 = fopen("Central_U_steady.plt","a+t");

        if ( fout2 == NULL )
        {
            printf("\nERROR when opening file\n");
            fclose( fout2 );
        }

        else
        {
            fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
            fprintf( fout2, "ZONE  T=\"Timestep 1\"\n");
            fprintf( fout2, "StrandID=1, SolutionTime=%d\n",step);
            fprintf( fout2, "I=%d, J=%d\n", grid, grid );

            for ( j = 0 ; j < (grid) ; j++ )
            {
                for ( i = 0 ; i < (grid) ; i++ )
                {
                    double xpos, ypos;
                    xpos = i*dx;
                    ypos = j*dy;

                    fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, u_c[i][j], v_c[i][j], p_c[i][j] );
                }
            }
        }

        fclose( fout2 );

        // CENTRAL --U
        fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
        fprintf( fout2, "ZONE  T=\"Timestep 1\"\n");
        fprintf( fout2, "StrandID=1, SolutionTime=%d\n",step);
        fprintf( fout2, "I=%d, J=%d\n");
        for ( j = 0 ; j < grid ; j++ )
        {
            double ypos;
            ypos = (double) j*dy;

            fprintf( fout3, "%5.8lf\t%5.8lf\n", (u_c[grid/2][j] + u_c[(grid/2)+1][j])/(2.), ypos  );
        }
        //step++;
     //}
     //return 0;
}

