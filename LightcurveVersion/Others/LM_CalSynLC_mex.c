#include "mex.h"
#include <math.h>
#include <stdio.h>
// SynLC = LM_CalSynLC_mex(NormalVec, Gw, E0, E, JDT, Q, RotatePole, G, cg9, cg10);
void CalSynLC(int TotalPoints, int Ndegree, double *NormalVec, double *Gw, double *E0, double *E, double *JDT, double *Q, double *RotatePole, double *G, double *cg9, double *cg10, double *SynLC)
{
    int i,f;
    double e0[3],e[3];
    double AsterE0[3], AsterE[3];
    double Temp[3],M[3];
    double len;
    double Phi;
    double cp,sp;
    double mu, mu0;
    
    for(i=0;i<TotalPoints;i++)
    {
        e0[0] = E0[i*3]; e0[1] = E0[i*3+1]; e0[2] = E0[i*3+2];
        e[0] = E[i*3];   e[1] = E[i*3+1];   e[2] = E[i*3+2];
        Phi = (*cg10) - (*cg9) *(JDT[i] - JDT[0]);
        cp = cos(Phi); sp = sin(Phi);
        // calculate Matrix * vectors: for E
        Temp[0] = RotatePole[0]*e[0] + RotatePole[3]*e[1] + RotatePole[6]*e[2];
        Temp[1] = RotatePole[1]*e[0] + RotatePole[4]*e[1] + RotatePole[7]*e[2];
        Temp[2] = RotatePole[2]*e[0] + RotatePole[5]*e[1] + RotatePole[8]*e[2];
        M[0] = cp * Temp[0] + sp * Temp[1];
        M[1] = -sp * Temp[0] + cp * Temp[1];
        M[2] = Temp[2];
        AsterE[0] = Q[0] * M[0] + Q[3] * M[1] + Q[6] * M[2] + G[0];
        AsterE[1] = Q[1] * M[0] + Q[4] * M[1] + Q[7] * M[2] + G[1];
        AsterE[2] = Q[2] * M[0] + Q[5] * M[1] + Q[8] * M[2] + G[2];
        // calculate Matrix * vectors: for E0
        Temp[0] = RotatePole[0]*e0[0] + RotatePole[3]*e0[1] + RotatePole[6]*e0[2];
        Temp[1] = RotatePole[1]*e0[0] + RotatePole[4]*e0[1] + RotatePole[7]*e0[2];
        Temp[2] = RotatePole[2]*e0[0] + RotatePole[5]*e0[1] + RotatePole[8]*e0[2];
        M[0] = cp * Temp[0] + sp * Temp[1];
        M[1] = -sp * Temp[0] + cp * Temp[1];
        M[2] = Temp[2];
        AsterE0[0] = Q[0] * M[0] + Q[3] * M[1] + Q[6] * M[2] + G[0];
        AsterE0[1] = Q[1] * M[0] + Q[4] * M[1] + Q[7] * M[2] + G[1];
        AsterE0[2] = Q[2] * M[0] + Q[5] * M[1] + Q[8] * M[2] + G[2];
        // normalizing
        len = sqrt(AsterE[0] * AsterE[0] + AsterE[1] * AsterE[1] + AsterE[2] * AsterE[2]);
        AsterE[0] /= len; AsterE[1] /= len; AsterE[2] /= len;
        
        len = sqrt(AsterE0[0] * AsterE0[0] + AsterE0[1] * AsterE0[1] + AsterE0[2] * AsterE0[2]);
        AsterE0[0] /= len; AsterE0[1] /= len; AsterE0[2] /= len;
        
        //  Calculate Total Brighness 
        for(f=0;f<Ndegree;f++)
        {
            mu = AsterE[0] * NormalVec[3*f] + AsterE[1] * NormalVec[3*f+1] + AsterE[2] * NormalVec[3*f+2];
            mu0 = AsterE0[0] * NormalVec[3*f] + AsterE0[1] * NormalVec[3*f+1] + AsterE0[2] * NormalVec[3*f+2];
            if(mu>0)
            {
                if(mu0>0)
                {
                    SynLC[i] += (mu*mu0*(0.1+1/(mu+mu0)))*Gw[f];
                }
            }
        }
    }
}
void mexFunction(int hlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *NormalVec, *Gw, *E, *E0, *JDT;
    double *Q, *RotatePole, *G, *cg9, *cg10;
    double *SynLC;
    int TotalPoints;
    int Ndegree;
    
    // fetch the pointers of input variables
    NormalVec = mxGetPr(prhs[0]);
    Gw = mxGetPr(prhs[1]);
    E0 = mxGetPr(prhs[2]);
    E = mxGetPr(prhs[3]);
    JDT = mxGetPr(prhs[4]);
    Q =  mxGetPr(prhs[5]);
    RotatePole = mxGetPr(prhs[6]);
    G =  mxGetPr(prhs[7]);
    cg9 =  mxGetPr(prhs[8]);
    cg10 =  mxGetPr(prhs[9]);
    
    TotalPoints = mxGetM(prhs[4]);
    Ndegree = mxGetM(prhs[1]);
    
    // create output variable
    plhs[0] = mxCreateDoubleMatrix(TotalPoints, 1, mxREAL);
    SynLC = mxGetPr(plhs[0]);
    // invoke the C function
    CalSynLC(TotalPoints, Ndegree,NormalVec, Gw, E0, E, JDT, Q, RotatePole, G, cg9, cg10, SynLC);
}
