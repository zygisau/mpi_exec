#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <array>

using namespace std;

int demandPointsCount = 1000;      // Vietoviu skaicius (demand points, max 10000)
int preexistingFacilitiesCount = 10;        // Esanciu objektu skaicius (preexisting facilities)
int candidateLocationsCount = 25;        // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 5;         // Nauju objektu skaicius

double **demandPoints; // Geografiniai duomenys


//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);
int increaseX(int *X, int index, int maxindex);
void splitProblemToParts(int* parts, int& numProcs, double& iterationTimes);
void computeAndCompare(int *X, int& numX, double& u, double& bestU, int *bestX);
double factorial(double n);

//=============================================================================


int main (int argc , char* argv[]) {
    MPI_Init(&argc, &argv);
    double ts = MPI_Wtime();

    loadDemandPoints();             // Nuskaitomi duomenys

    // Sudarom pradini sprendini: [0, 1, 2, 3, ...]
    int *X = new int[numX];
    int *bestX = new int[numX];
    for (int i=0; i<numX; i++) {
        X[i] = i;
        bestX[i] = i;
    }
    double u = evaluateSolution(X);
    double bestU = u;
    int r;
    //----- Pagrindinis ciklas ------------------------------------------------
    int MASTER_ID = 0;
    int id , numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Request request;
    MPI_Status mpiStatus;
    int* parts = new int[numProcs];
    double iterationTimes = factorial(candidateLocationsCount) / (factorial(numX) * factorial(candidateLocationsCount - numX));
    splitProblemToParts(parts, numProcs, iterationTimes);
    bool stopAlgo;
    if (id == MASTER_ID) {
        for (int i = 0; i < iterationTimes; ++i) {
            for (int j = 1; j < numProcs; ++j) {
                if (increaseX(X, numX - 1, candidateLocationsCount)) {
                    int tag = j;
                    MPI_Isend(X, numX, MPI_INT, j, tag, MPI_COMM_WORLD, &request);
                    i++;
                } else {
                    stopAlgo = true;
                    break;
                }
            }
            if (stopAlgo) {
                break;
            }

            if (increaseX(X, numX - 1, candidateLocationsCount)) {
                computeAndCompare(X, numX, u, bestU, bestX);
            } else {
                stopAlgo = true;
                break;
            }
        }
    } else {
        for (int i = 0; i < parts[id]; ++i) {
            int tag = id;
            MPI_Recv(X, numX, MPI_INT, MASTER_ID, tag, MPI_COMM_WORLD, &mpiStatus);
            computeAndCompare(X, numX, u, bestU, bestX);
        }
    }

    if (id == MASTER_ID) {
        int sentU;
        int *sentX = new int[numX];
        for (int i = 1; i < numProcs; ++i) {
            MPI_Recv(&sentU, 1, MPI_INT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &mpiStatus);
            int tag = mpiStatus.MPI_SOURCE;
            if (sentU > bestU) {
                bestU = sentU;
                MPI_Recv(sentX, numX, MPI_INT, MPI_ANY_SOURCE, tag*10, MPI_COMM_WORLD, &mpiStatus);
                for (int i=0; i<numX; i++) bestX[i] = sentX[i];
            }
        }
    } else {
        MPI_Send(&bestU, 1, MPI_INT, MASTER_ID, id, MPI_COMM_WORLD);
        MPI_Isend(bestX, numX, MPI_INT, MASTER_ID, id*10, MPI_COMM_WORLD, &request);
    }

    //----- Rezultatu spausdinimas --------------------------------------------
    if (id == MASTER_ID) {
        cout << "Geriausias sprendinys: ";
        for (int i=0; i<numX; i++) cout << bestX[i] << " ";
        cout << "(" << bestU << ")" << endl;
    }
    double tf = MPI_Wtime();
    cout << "Procesorius #" << id << ";\t" << "Skaiciavimo trukme: " << tf-ts << " s." << endl;

    MPI_Finalize();
}

void computeAndCompare(int *X, int& numX, double& u, double& bestU, int *bestX) {
    u = evaluateSolution(X);
    if (u > bestU) {
        bestU = u;
        for (int i=0; i<numX; i++) bestX[i] = X[i];
    }
}

double factorial(double n) {
    if (n < 0) {
        return 0;
    }
    return !n ? 1.0 : (n * factorial(n - 1));
}

void splitProblemToParts(int* parts, int& numProcs, double& iterationTimes) {
    for (int i = 0; i < numProcs; ++i) {
        parts[i] = iterationTimes / numProcs;
    }
    for (int i = 0; i < candidateLocationsCount%numProcs; ++i) {
        parts[i] += candidateLocationsCount%numProcs;
    }
}

//=============================================================================

void loadDemandPoints() {
    FILE *f;
    f = fopen("demandPoints.dat", "r");
    demandPoints = new double*[demandPointsCount];
    for (int i=0; i < demandPointsCount; i++) {
        demandPoints[i] = new double[3];
        fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
    }
    fclose(f);
}

//=============================================================================

double HaversineDistance(double* a, double* b) {
    double dlon = fabs(a[0] - b[0]);
    double dlat = fabs(a[1] - b[1]);
    double aa = pow((sin((double)dlon/(double)2*0.01745)),2) + cos(a[0]*0.01745) * cos(b[0]*0.01745) * pow((sin((double)dlat/(double)2*0.01745)),2);
    double c = 2 * atan2(sqrt(aa), sqrt(1-aa));
    double d = 6371 * c;
    return d;
}

//=============================================================================

double getTime() {
    struct timeval laikas;
    gettimeofday(&laikas, NULL);
    double rez = (double)laikas.tv_sec+(double)laikas.tv_usec/1000000;
    return rez;
}

//=============================================================================

double evaluateSolution(int *X) {
    double U = 0;
    int bestPF;
    int bestX;
    double d;
    for (int i=0; i < demandPointsCount; i++) {
        bestPF = 1e5;
        for (int j=0; j < preexistingFacilitiesCount; j++) {
            d = HaversineDistance(demandPoints[i], demandPoints[j]);
            if (d < bestPF) bestPF = d;
        }
        bestX = 1e5;
        for (int j=0; j<numX; j++) {
            d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
            if (d < bestX) bestX = d;
        }
        if (bestX < bestPF) U += demandPoints[i][2];
        else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
    }
    return U;
}

//=============================================================================

bool isNextLocationAvailable(int *X, int index, int maxCityIndex) {
    int nextLocationIndex = X[index]+1;
    int occupiedIndexCount = numX - index - 1;
    int availableNewCitiesMaxIndex = maxCityIndex - occupiedIndexCount;
    return nextLocationIndex < availableNewCitiesMaxIndex;
}

int increaseX(int *X, int index, int maxindex) {
//    if next location is available
    if (isNextLocationAvailable(X, index, maxindex)) {
        X[index]++;
    }
    else {
        bool isLastCity = index == 0;
        bool isNextLocationMaxAvailableLocationsIndex = X[index]+1 == maxindex-(numX-index-1);
        if (isLastCity && isNextLocationMaxAvailableLocationsIndex) {
            return 0;
        }
        else {
            if (increaseX(X, index-1, maxindex)) X[index] = X[index-1]+1;
            else return 0;
        }
    }
    return 1;
}