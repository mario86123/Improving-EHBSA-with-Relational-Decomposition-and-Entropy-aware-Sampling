
#include "VRP.h"
#include "Tools.h"

// #define for_(i, a, b) for (int i=(a);i<(b);++i)
#define sqr(a) (a)*(a)

/*
 *Class constructor.
 */
VRP::VRP()
{
	
}

/*
 * Class destructor.
 */
VRP::~VRP()
{
    for (int i = 0; i < m_size; i++)
        delete [] cost[i];
    delete [] cost;
    delete [] demand;

}



int VRP::Read(string filename)
{
    string line;
    ifstream file(filename);
        
    while (file) {
        getline(file, line);
        
        if (line.empty()) continue;
        
        string key_word, x;
        istringstream ss(line);
        
        ss >> key_word;
        
        if (key_word == "DIMENSION") {
            
            ss >> x >> numCustomer;
            
            cost = new double * [numCustomer];
            for (int i = 0; i < numCustomer; i++)
            {
                cost[i]= new double[numCustomer];
            }

            demand = new int[numCustomer];
        }

        else if (key_word == "CAPACITY") {
            ss >> x >> capacity;
        }

        else if (key_word == "TRUCKS") {
            ss >> x >> numTrucks;
            m_size = numTrucks + numCustomer - 2;
            // cout << "numCustomer: " << numCustomer << endl; 
            // cout << "numTrucks: " << numTrucks << endl; 

        }

        else if (key_word == "NODE_COORD_SECTION") {
            
            vector<int> x(numCustomer), y(numCustomer);
            
            int a, b, c;
            for (;;) {
                getline(file, line);
                
                if (line[0] == 'D')
                break;
                
                istringstream ss(line);
                
                ss >> a >> b >> c;
                
                x[a-1] = b;
                y[a-1] = c;
            }
            
            for(int i = 0; i < numCustomer; i++) {
                for(int j = 0; j < numCustomer; j++) {
                    double dist2 = sqr(x[i]-x[j]) + sqr(y[i]-y[j]);
                    double dist = sqrt(dist2);
                    // cout << "i: " << i << ", j: " << j << " ==> " << dist << endl;
                    cost[i][j] = floor(dist + 0.5);
                    cost[j][i] = floor(dist + 0.5);
                }
            }
            
            // PrintMatrix(cost, numCustomer);
            // DEMAND_SECTION
            
            for (;;) {
                getline(file, line);
                
                if (line[0] == 'D')
                break;
                
                istringstream ss(line);
                int v, d;
                ss >> v >> d;
                demand[v-1] = d;
            }
        
        }
    }
    // cout << "m_size: " << m_size << endl; 
	return m_size;
}


#define K 1
/*
 * This function evaluates the individuals for the QAP problem.
 */
double VRP::Evaluate(int * genes)
{
    stack<int> route;  
    int* weight = new int[numTrucks];
    double* length = new double[numTrucks];
    
    int vehicle = 0;
    
    while(!route.empty())
        route.pop();

    for(int i = 0; i <= m_size; i++)
    {
        int index;

        // add genes[i] to route stack
        if (i < m_size){
            index = genes[i] + 1;
        }
        else { // add terminal node to route stack
            index = m_size + 1;
        }
        // cout << "index: " << index << endl;

        if (index < numCustomer) {
            // route.push(index);
            route.push(index);
        }
        else 
        {
            // cout << "index: " << index << endl;
            if (route.empty())
            {
                length[vehicle] = 0.0;
                weight[vehicle] = 0;
            } // if
            else
            {
                int w = 0;
                double l = 0.0;
                int top = 0;

                while(!route.empty())
                {
                    int next = route.top();
                    route.pop();
                    // cout << "next: " << next << ", top: " << top << ", cost: " << cost[next][top] << endl;
                    l += cost[next][top];                 
                    w += demand[next]; 
                    top = next;            
                } // while
                l += cost[0][top];                
                
                // printf("l: %lf\n", l);
                length[vehicle] = l;
                weight[vehicle] = w;

                vehicle++;     
                // PrintArray(length, numTrucks);           
                // PrintArray(weight, numTrucks);           
            } // else
        } // else
    } // for
    
    double fitness = 0.0;
    for (int i = 0; i < vehicle; i++)
    {
        double exceed = (double)(weight[i] - capacity) / ((double)capacity);
        double penalty = K * length[i] * exceed;
        double fitness_i = length[i] + max(0.0, penalty);
        // cout << " max(0.0, penalty): " << max(0.0, penalty) << endl;
        fitness += fitness_i;
        // cout << i << ") length[i]: "<< length[i] <<", penalty: " << penalty << ", fitness_i: " << fitness_i << endl;
    }
    // cout << endl;
    return -fitness;

}

/*
 * This function evaluates the inverted solution of the given individual for the QAP problem.
 */
double VRP::EvaluateInv(int * genes)
{
    return 0;
}

/*
 * Returns the size of the problem.
 */
int VRP::GetProblemSize()
{
    return m_size;
}