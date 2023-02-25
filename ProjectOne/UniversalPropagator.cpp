//Universal Variable two body orbit propagator, Curtis ALgorithm 3.3 and 3.4

#include <vector>
#include <cmath>
#include <iostream>
#include "UniversalPropagator.h"
using namespace std; // Don't need to use std:: everywhere 

const double pi = acos(-1);


double Stumpff_C(const double z) {
    if (z > 0) {
        return (1-cos(sqrt(z)))/z;
    } else if (z < 0) {
        return (cosh(sqrt(-z))-1)/(-z);
    } else {
        return 1/2.0;
    }
}

double Stumpff_S(const double z) {
    if (z > 0) {
        return (sqrt(z)-sin(sqrt(z)))/pow(sqrt(z),3);
    } else if (z < 0) {
        return (sinh(sqrt(-z))-sqrt(-z))/pow(sqrt(-z),3);
    } else {
        return 1/6.0;
    }
}

double solveUniversalAnomaly(double r0, double v0, double mu, double alpha, double dt) {
    double chi =sqrt(mu)*alpha*dt;
    double chiNew = chi;
    double tol = 1e-8;
    double ratio = 1;
    double z = 0; double f = 0; double fPrime = 0;  ;
    int iter = 0;
    while ((ratio > tol) && (iter < 100)) {
        z = alpha*pow(chi,2);
        f = r0*v0*pow(chi,2)*Stumpff_C(z)/sqrt(mu) + (1-alpha*r0)*pow(chi,3)*Stumpff_S(z) + r0*chi - sqrt(mu)*dt;
        fPrime = r0*v0*chi*(1-z*Stumpff_S(z))/sqrt(mu) + (1-alpha*r0)*pow(chi,2)*Stumpff_C(z) + r0;
        ratio = fabs(f/fPrime);
        chiNew = chi - f/fPrime;
        chi = chiNew;        
        iter++;
    }
    if (iter==100) {
        cout<<"Universal Propagator did not converge in 100 iterations"<<endl;
    } 
    return chiNew;
}

double norm(const std::vector<double>& v) {
    double norm = 0;
    if (v.size() == 0) {
        cout << "Vector is empty" << endl;
    }
    for (int i=0; i<v.size(); i++) {
        norm += pow(v[i],2);
    }
    return sqrt(norm);
}

std::vector<double> multiply_vector(const std::vector<double>& v, const double c) {
    if (v.size() == 0) {
        cout << "Vector is empty" << endl;
    }
    std::vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++) {
        result[i] = v[i] * c;
    }
    return result;
}

std::vector<double> add_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result(v1.size());
    if (v1.size() != v2.size()) {
        cout << "Vectors are not the same size in add_vectors" << endl;
    }
    if (v1.size() == 0) {
        cout << "Vector is empty" << endl;
    }
    for (int i = 0; i < v1.size(); i++) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

std::vector <double> subtract_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result(v1.size());
    if (v1.size() != v2.size()) {
        cout << "Vectors are not the same size in subtract_vectors" << endl;
    }
    if (v1.size() == 0) {
        cout << "Vector is empty" << endl;
    }
    for (int i = 0; i < v1.size(); i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

std::vector<double> CrossProduct(std::vector<double>& v1, std::vector<double>& v2) {
    std::vector<double> result(3);
    if (v1.size() != 3 || v2.size() != 3) {
        cout << "Vectors are not 3D" << endl;
    }
    result[0] = v1[1]*v2[2] - v1[2]*v2[1];
    result[1] = v1[2]*v2[0] - v1[0]*v2[2];
    result[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return result;
}

double dotProduct(std::vector<double>& v1, std::vector<double>& v2) {
    double result = 0;
    if (v1.size() != v2.size()) {
        cout << "Vectors are not the same size" << endl;
    }
    if (v1.size() == 0) {
        cout << "Vector is empty" << endl;
    }
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i]*v2[i];
    }
    return result;
}

double minValueInArray(std::vector<double> array) {
    if (array.size() == 0) {
        return 0.0;
        cout << "Array is empty" << endl;
    }
    else {
        double min = array[0];
        for (int i = 1; i < array.size(); i++) {
            if (array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }
} 

std::vector<double> GetPositionVelocity(std::vector<double> r0, std::vector<double> v0, double mu, double dt) {
    double norm_r0 = norm(r0);
    double norm_v0 = norm(v0);
    double projvr0 = dotProduct(r0,v0)/norm_r0; // projection of v0 onto r0 (scalar
    double alpha = 2/norm_r0 - pow(norm_v0,2)/mu;
    double chi = solveUniversalAnomaly(norm_r0, projvr0, mu, alpha, dt);
    double z = alpha*pow(chi,2);
    double f = 1-pow(chi,2)/norm_r0*Stumpff_C(z);
    double g = dt-pow(chi,3)/sqrt(mu)*Stumpff_S(z);
    
    std::vector<double> r = add_vectors(multiply_vector(r0,f), multiply_vector(v0,g));
    
    double norm_r = norm(r);
    double fDot = sqrt(mu)*chi*(z*Stumpff_S(z)-1)/(norm_r*norm_r0);
    double gDot = 1-pow(chi,2)*Stumpff_C(z)/norm_r;

    std::vector<double> v = add_vectors(multiply_vector(r0,fDot), multiply_vector(v0,gDot));
    std::vector<double> rv(6, 0.0);
    rv[0]= r[0];
    rv[1]= r[1];
    rv[2]= r[2];
    rv[3]= v[0];
    rv[4]= v[1];
    rv[5]= v[2];

    return rv;
}

double hypergeometricF(double z, double tol)
// Code was taken from Izzo's Github Repo
{
    double Sj = 1.0;
    double Cj = 1.0;
    double err = 1.0;
    double Cj1 = 0.0;
    double Sj1 = 0.0;
    int j = 0;
    while (err > tol) {
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
        Sj1 = Sj + Cj1;
        err = fabs(Cj1);
        Sj = Sj1;
        Cj = Cj1;
        j = j + 1;
    }
    return Sj;
}

double timeOfFlight(const double x, const int M, const double lambda) {
    double dist = fabs(x-1);
    double a_amin = 1.0/(1-x*x); // ratio a/amin
    //cout << "a_amin: " << a_amin << endl;
    //cout << "dist: " << dist << endl;
    double T_nonDim = 0;
    
    double battin = 0.01;
    double lagrange = 0.2;
    if (dist < lagrange && dist > battin) { 
        if (a_amin>0) {
            double alpha = 2.0*acos(x); // eq 14
            double beta = 2.0 * asin(sqrt(lambda*lambda/a_amin)); // eq 13
            //TOF equation 9, but non dimensionalized (see eq for T under eq 13)
            T_nonDim = a_amin*sqrt(a_amin)*(alpha - beta + sin(beta) - sin(alpha) + 2*M*pi)/2;
        } else {
            double alpha = 2.0 * acosh(x);
            double beta = 2.0 * asinh(sinh(alpha/2.0)/lambda); 
            T_nonDim = (-a_amin * sqrt(-a_amin) * ((beta - sinh(beta)) - (alpha - sinh(alpha))) / 2.0);; 

        }
        return T_nonDim;
    }
    // Code taken from Izzo Github Repo: https://github.com/esa/pykep/blob/master/src/lambert_problem.cpp
    double K = lambda * lambda;
    double E = x * x - 1.0;
    double rho = fabs(E);
    double z = sqrt(1 + K * E);
    if (dist < battin) { // We use Battin series tof expression
        double eta = z - lambda * x;
        double S1 = 0.5 * (1.0 - lambda - x * eta);
        double Q = hypergeometricF(S1, 1e-11);
        Q = 4.0 / 3.0 * Q;
        T_nonDim = (eta * eta * eta * Q + 4.0 * lambda * eta) / 2.0 + M * M_PI / pow(rho, 1.5);
        return T_nonDim;
    } else { // We use Lancaster tof expresion
        double y = sqrt(rho);
        double g = x * z - lambda * E;
        double d = 0.0;
        if (E < 0) {
            double l = acos(g);
            d = M * M_PI + l;
        } else {
            double f = y * (z - lambda * x);
            d = log(f + g);
        }
        T_nonDim = (x - lambda * z - d / y) / E;
        return T_nonDim;
    }
}

void dTdx(double &DT, double &DDT, double &DDDT, const double x, const double T, const double lambda){
    double l2 = lambda * lambda;
    double l3 = l2 * lambda;
    double umx2 = 1.0 - x * x;
    double y = sqrt(1.0 - l2 * umx2);
    double y2 = y * y;
    double y3 = y2 * y;
    DT = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y);
    DDT = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3);
    DDDT = 1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / y3 / y2);
}

int HouseholderIterations(double &x0, const double T, const int N, const double tol, const int maxIter, const double lambda) {
    int iter = 0; 
    double error = 1; 
    double DT = 0; double DDT = 0; double DDDT = 0;
    double xNew = 0; double f = 0; double Tout = 0; //Initialize
    
    while ((iter <= maxIter) && (error > tol)) {
        Tout = timeOfFlight(x0, N, lambda);
        dTdx(DT, DDT, DDDT, x0, T, lambda);
        f = Tout - T; 
        xNew = x0 - f * (pow(DT,2) - f * DDT / 2.0) / (DT * (pow(DT,2) - f * DDT) + DDDT*f*f/6.0);                  
        error = fabs(xNew - x0);            
        x0 = xNew;
        iter++;
    }

    if (iter > maxIter) {
        cout << "Householder iterations did not converge, x = " << x0 << endl;
    }
    return iter;
}

std::vector<double> Get_x_ListIzzo(double lambda, double T){
    int m_Max = floor(T/pi);
    double T0 = acos(lambda) + lambda*sqrt(1-pow(lambda,2)); // single revolution case 
    double TM0 = T0 + m_Max * M_PI; 
    double T1 = 2.0*(1-pow(lambda,3))/3.0;
    if (m_Max > 0) { // Do halley iterations to determine if there is
    // enough time for an elliptic arc with M revolutions that hits the final point
        if (T < TM0) {
            int iter = 0; 
            double error = 1; 
            double T_min = TM0; 
            double x_old = 0; double x_new = 0; 
            double DT = 0; double DDT = 0; double DDDT = 0;
            while (1) {
                // compute derivatives, eq22 of Izzo paper. 
                dTdx(DT, DDT, DDDT, x_old, T, lambda);
                if (DT !=0 ){
                    x_new = x_old - (DT*DDT)/(DDT*DDT - 0.5*DT*DDDT);
                }
                error = fabs(x_new - x_old);

                if ((error < 1e-10) || (iter == 15)) {
                    break;
                    // Compute time of flight from x_new
                    T_min = timeOfFlight(x_new, m_Max, lambda);
                }   
                x_old = x_new;
                iter++;
            }
            if (T_min > T) {
                m_Max = m_Max - 1;
            } 
        }
    }
    // Initialize vector output for x and y 
    // Typically 4 ellipses for each Lambert problem, 2 prograde and 2 retrograde. Here we only do one direction (prograde)
    std::vector<int> x_iter(m_Max*2+1); 
    std::vector<double> x(m_Max*2+1); 
    std::vector<double> y(m_Max*2+1);
    double tol = 1e-8; 
    // Single revolution Case - initial guess, eq 30
    if (T >= T0) { 
        x[0] = pow((T0/T),2/3)-1;
    } else if (T <= T1) {
        x[0] = 2*(T1/T) - 1; //2.5 * T1 * (T1 - T) / ((1 - pow(lambda,5)) * T) + 1;
    } else {
        x[0] = pow((T/T0), 0.69314718055994529 / log(T1/T0)) - 1.0;
    }
    bool debug = false; 
    if (debug) {
        cout << "T = " << T << endl;
        cout << "T0 = " << T0 << endl;
        cout << "T1 = " << T1 << endl;
        cout << "TM0 = " << TM0 << endl;
        cout << "lambda = " << lambda << endl;
        cout << "x[0] = " << x[0] << endl;
    } 
    // Householder iterations to find x
    x_iter[0] = HouseholderIterations(x[0], T, 0, tol, 25, lambda);
    //cout << "x[0] = " << x[0] << endl;
    //Multiple revolution case 
    double temp; 
    for (int M = 1; M <= m_Max; M++) { 
        temp = pow((M * M_PI + M_PI) / (8.0 * T), 2.0 / 3.0);
            //cout << "Revoultion Case = " << M << endl;
            // 3.2.1 left Householder iterations
            x[2 * M - 1] = (temp - 1) / (temp + 1);
            x_iter[2 * M - 1] = HouseholderIterations(x[2 * M - 1], T, M, tol, 15, lambda);
            // 3.2.1 right Householder iterations
            temp = pow((8.0 * T) / (M * M_PI), 2.0 / 3.0);
            x[2 * M] = (temp - 1) / (temp + 1);
            x_iter[2 * M] = HouseholderIterations(x[2 * M], T, M, tol, 15, lambda);
    }
    return x;

} // end of findxListIzzo

std::vector<double> LambertSolverIzzo(std::vector<double>& rv1, std::vector<double>& rv2, double mu, double timeOfFlight, const bool cw, const bool rendezvous){
    // Izzo Lambert Solver
    std::vector <double> r1 = {rv1[0], rv1[1], rv1[2]};
    std::vector <double> r2 = {rv2[0], rv2[1], rv2[2]};
    std::vector <double> vDepBody = {rv1[3], rv1[4], rv1[5]};
    std::vector <double> vArrBody = {rv2[3], rv2[4], rv2[5]};
    
    double r1_norm = norm(r1);
    double r2_norm = norm(r2);
 
    double c =  norm(subtract_vectors(r2,r1)); 
    double s = (r1_norm+r2_norm+c)/2;
    double lambda = sqrt(1-c/s);
    double T = sqrt(2.0 * mu /pow(s,3)) * timeOfFlight;
    
    std::vector<double> ir1 = multiply_vector(r1,1/r1_norm);
    std::vector<double> ir2 = multiply_vector(r2,1/r2_norm);
    std::vector<double> ih = CrossProduct(ir1,ir2);
    std::vector<double> it1;
    std::vector<double> it2;

    ih = multiply_vector(ih,1/norm(ih));
    
    if (ih[2] == 0) {
        cout << "Angular momentum has no z component" << endl;
    }
    if (ih[2] < 0) { // Transfer angle is larger than 180 degrees
        lambda = -lambda;
        it1 = CrossProduct(ir1,ih);
        it2 = CrossProduct(ir2,ih);
    } else {
        it1 = CrossProduct(ih,ir1);
        it2 = CrossProduct(ih,ir2);
    }

    it1 = multiply_vector(it1,1/norm(it1));
    it2 = multiply_vector(it2,1/norm(it2));

    if (cw) { // Retrograde motion
        lambda = -lambda;
        it1[0] = -it1[0];
        it1[1] = -it1[1];
        it1[2] = -it1[2];
        it2[0] = -it2[0];
        it2[1] = -it2[1];
        it2[2] = -it2[2];
    }

    // Calculate maximum possible value of M (number of revolutions) using Halley's Algorithm
    // Get list of x,y values for each M using Householder iterations
    std::vector<double> x = Get_x_ListIzzo(lambda, T);
    int numSolutions = x.size();
    
    double gamma = sqrt(mu*s/2); 
    double rho = (r1_norm-r2_norm)/c;
    double sigma = sqrt(1-pow(rho,2));
    double y = 0; 
    double l2 = lambda * lambda;
    double Vr1 = 0; double Vr2 = 0; double Vt1 = 0; double Vt2 = 0;
    std::vector <double> V1; std::vector <double> V2;
    std::vector<double> totalDV(numSolutions, 0.0);

    for (int i = 0; i < numSolutions; i++) {
        // Calculate delta velocity vectors for each x value
        y = sqrt(1-l2 + l2*x[i]*x[i]);
        
        Vr1 = gamma * (lambda*y-x[i]-rho*(lambda*y+x[i]))/r1_norm;
        Vr2 = -gamma * (lambda*y-x[i]+rho*(lambda*y+x[i]))/r2_norm;
        Vt1 = gamma * sigma * (y +lambda*x[i])/r1_norm;
        Vt2 = gamma * sigma * (y +lambda*x[i])/r2_norm;

        V1 = add_vectors(multiply_vector(ir1,Vr1),multiply_vector(it1,Vt1));
        V2 = add_vectors(multiply_vector(ir2,Vr2),multiply_vector(it2,Vt2));
        if (rendezvous) {
            totalDV[i] = norm(subtract_vectors(V1,vDepBody)) + norm(subtract_vectors(V2,vArrBody));
        } else{
            totalDV[i] = norm(subtract_vectors(V1,vDepBody));
        }
    }
    std::vector<double> out = {V1[0], V1[1], V1[2], V2[0], V2[1], V2[2]};
    //return out;
    return totalDV;
}

std::vector<std::vector<double>> AssignmentOne(int targetNum, bool rendezvous) {
    // Problem Definition 
    std::vector<double> r1I = {3.515868886595499*pow(10,-2), -3.162046390773074, 4.493983111703389};
    std::vector<double> v1I = {-2.317577766980901*pow(10,-3), 9.843360903693031*pow(10,-3), -1.541856855538041*pow(10,-2)};
    std::vector<double> r2I = {7.249472033259724, 14.61063037906177, 14.24274452216359};
    std::vector<double> v2I = {-8.241709369476881*pow(10,-3), -1.156219024581502 *pow(10,-2), -1.317135977481448 *pow(10,-2)};
    std::vector<double> rE = {-1.796136509111975*pow(10,-1), 9.667949206859814 *pow(10,-1), -3.668681017942158 *pow(10,-5)};
    std::vector<double> vE = {-1.720038360888334*pow(10,-2), -3.211186197806460 *pow(10,-3), 7.927736735960840 *pow(10,-7)};
    
    //Constants 
    double AU = 149597870.7; //km
    double day = 86400; //seconds
    //mu of sun in units AU, au/day
    double mu_sun = 1.32712440018*pow(10,11);
    mu_sun = mu_sun/pow(AU,3) * pow(day,2);
    
    //Initialize Outputs
    std::vector<double> rvEarthDepartureTime(6,0.0);
    std::vector<double> rvTargetAtArrival(6,0.0);
    std::vector<double> targetPosDeparture(3,0.0);
    std::vector<double> targetVelDeparture(3,0.0);
    
    double departureWindowLength = 365; //days
    double arrivalWindowLength = 182+365; //days Aug 2017 - Jan 2019
    double initialDepartureDay = 1; 
    double initialArrivalDay = 185;

    if (targetNum == 1) {
        targetPosDeparture = r1I;
        targetVelDeparture = v1I;
    } else {
        targetPosDeparture = r2I;
        targetVelDeparture = v2I;
        departureWindowLength = 365*3+182; // = days between Jan2017 - July 2020
        arrivalWindowLength = 210+365*2; // = days between June 2019 - Jan 2022
        initialDepartureDay = 1; // Days after ephemerides (Jan 2017)
        initialArrivalDay = 185+365*2; // June 2019 - Jan 2017
    }   

    // cast double to int
    int numDepartureTimes = (int)departureWindowLength + 1;
    int numArrivalTimes = (int)arrivalWindowLength + 1;
    
    std::vector<double> departureTimes(numDepartureTimes,0.0);
    std::vector<double> arrivalTimes(numArrivalTimes,0.0);
    std::vector<std::vector<double>> totalDV(numDepartureTimes, std::vector<double>(numArrivalTimes, 0.0));

    bool useSI = true;
    if (useSI) {
        targetPosDeparture = multiply_vector(targetPosDeparture,AU);
        targetVelDeparture = multiply_vector(targetVelDeparture,AU/day);
        rE = multiply_vector(rE,AU);
        vE = multiply_vector(vE,AU/day);
        departureWindowLength = departureWindowLength*day;
        arrivalWindowLength = arrivalWindowLength*day;
        initialDepartureDay = initialDepartureDay*day;
        initialArrivalDay = initialArrivalDay*day;
        mu_sun = 1.32712440018e11; //km^3/s^2
    } 

    for (int i = 0; i < numDepartureTimes; i++) {
        departureTimes[i] = initialDepartureDay + i * departureWindowLength/(numDepartureTimes-1);
        //cout << "Departure Time: " << departureTimes[i] << endl;
    }
    for (int i = 0; i < numArrivalTimes; i++) {
        arrivalTimes[i] = initialArrivalDay + i * arrivalWindowLength/(numArrivalTimes-1);
        //cout << "Arrival Time: " << arrivalTimes[i] << endl;
    }

    double dt = 0;
    
    for (int i = 0; i < numDepartureTimes; i++) {
        for (int j = 0; j < numArrivalTimes; j++) {
            dt = arrivalTimes[j] - departureTimes[i];
            if (dt>0) {
                rvEarthDepartureTime = GetPositionVelocity(rE,vE,mu_sun,departureTimes[i]);
                rvTargetAtArrival = GetPositionVelocity(targetPosDeparture,targetVelDeparture,mu_sun,arrivalTimes[j]);

                std::vector<double> possibleDV = LambertSolverIzzo(rvEarthDepartureTime, rvTargetAtArrival, mu_sun, dt, false, rendezvous);
                totalDV[i][j] = minValueInArray(possibleDV);                
            }
            else {
                totalDV[i][j] = NAN;
            }
        }
    }
    return totalDV;
}

int main() {

    //Testing propagator and Lambert Solver 

    std::vector<double> r1 = {5644, 2830, 4170, 0, 0, 0};
    std::vector<double> r2 = {-2240, 7320, 4980, 0, 0, 0};
    double mu = 3.986004418e5;
    double tof = 1200;

    bool cw = false; 
    std::vector<double> out = LambertSolverIzzo(r1, r2, mu,  tof, cw, true);
    cout << "out: " << out[0] << ", " << out[1] << ", " << out[2] << endl;
    cout << "out: " << out[3] << ", " << out[4] << ", " << out[5] << endl;

    cw = true; 
    out = LambertSolverIzzo(r1, r2, mu,  tof, cw, true);
    cout << "out: " << out[0] << ", " << out[1] << ", " << out[2] << endl;
    cout << "out: " << out[3] << ", " << out[4] << ", " << out[5] << endl;

    std::vector<double> v1 = {out[0], out[1], out[2]};
    std::vector<double> r1Check = {r1[0], r1[1], r1[2]};
    std::vector<double> rvCheck = GetPositionVelocity(r1Check, v1, mu, tof);
    //cout << "posFinal: " << rvCheck[0] << ", " << rvCheck[1] << ", " << rvCheck[2] << endl;
    //cout << "velFinal: " << rvCheck[3] << ", " << rvCheck[4] << ", " << rvCheck[5] << endl;

    // Check Lambert problem on Slack: 
    std::vector<double> rRetro = {5644, 2830, 4170};
    std::vector<double> vRetro = {-6.11037079, -5.96890258, -6.79748569};
    std::vector<double> rvCheckRetrograde2 = GetPositionVelocity(rRetro, vRetro, mu, tof);
    //cout << "posFinal: " << rvCheckRetrograde2[0] << ", " << rvCheckRetrograde2[1] << ", " << rvCheckRetrograde2[2] << endl;
    //cout << "velFinal: " << rvCheckRetrograde2[3] << ", " << rvCheckRetrograde2[4] << ", " << rvCheckRetrograde2[5] << endl;

    //Test Universal propagator - Example 3.7 from Curtis
    r1Check = {7000.0,-12124.0,0.0};
    v1 = {2.6679,4.6210,0.0};
    //std::vector<double> rvCheck2 = GetPositionVelocity(r1Check, v1, mu, 3600.0);
    //cout << "posFinal: " << rvCheck2[0] << ", " << rvCheck2[1] << ", " << rvCheck2[2] << endl;
    //cout << "velFinal: " << rvCheck2[3] << ", " << rvCheck2[4] << ", " << rvCheck2[5] << endl;

    //double chi = solveUniversalAnomaly(10000, 3.0752, 398600, -5.0878e-5, 3600);
    //cout << "chi: " << chi << endl;
}