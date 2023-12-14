#include <omp.h>
#include <chrono> 
#include <iostream> 
#include <string>
#include <vector>
#include <immintrin.h>
#include <sstream>
#include <fstream>
#include "quantum.cpp"
#include "gates.cpp"
#include <bitset>
//#include "tableau.cpp"

#define         CNOT		0
#define         HADAMARD	1
#define         PHASE		2
#define         MEASURE		3


void print_m512i(__m512i var) {
    uint64_t val[8];
    _mm512_storeu_si512(val, var);
    for (int i = 0; i < 8; i++) {
        std::bitset<64> bits(val[i]);
        std::cout << bits << std::endl;
    }
    std::cout << "\n";
}
void runprog(struct QProg *h, struct QState *q)

// Simulate the quantum circuit

{

	int m; // measurement result
	char not_measured = 1;

	for (int t = 0; t < h->gate_count; t++)
	{
         //if (h->op[t]==CNOT) cnot(q,h->control_qubit[t],h->target_qubit[t]);
         if (h->op[t]==HADAMARD) hadamard(q,h->control_qubit[t]);
         if (h->op[t]==PHASE) phase(q,h->control_qubit[t]);
        //  if (h->op[t]==MEASURE)
        //  {
        //          not_measured = 0;
        //          m = measure(q,h->control_qubit[t],h->SUPPRESSM);
        //          if (!h->SILENT)
        //          {
        //                  printf("\nOutcome of measuring qubit %ld: ", h->b[t]);
        //                  if (m>1) printf("%d (random)", m-2);
        //                  else printf("%d", m);
        //          }
        //  }
	}
	printf("\n");
	// if (h->DISPTIME)
	// {
    //      dt = difftime(time(0),tp);
    //      printf("\nMeasurement time: %lf seconds", dt);
    //      printf("\nTime per 10000 measurements: %lf seconds\n", dt*10000.0f/h->num_qubits);
	// }
	return;

}

/**
 * Initializes the state of a quantum system.
 * 
 * @param q A pointer to the QState structure representing the quantum system.
 * @param n The number of qubits in the system.
 * @param s A character array representing the state of the system.
 */
void initialize_state(QState *q, int n, char *s)
{
    
    q->num_qubits = n;
    q->x_bits = new __m512i*[(2*n + 1)];
    q->z_bits = new __m512i*[(2*n + 1)];
    q->phase = new int[2*n + 1];
    q->pw[0] = 1;
    q->over512 = (n/512)+1; // Changed from (n/32)+1 to (n/512)+1
    for (int i = 1; i < 32; i++)
        q->pw[i] = 2*q->pw[i-1];
    for (int i = 0; i < 2*n + 1; i++)
    {
        q->x_bits[i] = new __m512i[q->over512]; // Changed from over32 to over512
        q->z_bits[i] = new __m512i[q->over512]; // Changed from over32 to over512
        for (int j = 0; j < q->over512; j++)
        {
            //std::cout << "i: " << i << " j: " << j << "\n";
            q->x_bits[i][j] = _mm512_set1_epi64(0);
            q->z_bits[i][j] = _mm512_set1_epi64(0);
        }
        if (i < n){
            q->x_bits[i][i>>6] = _mm512_set1_epi64(q->pw[i&63]); // Changed from i>>5 to i>>6 and i&31 to i&63
            //print_m512i(q->x_bits[i][i>>6]);
        }
        else if (i < 2*n) {
            q->z_bits[i][(i-n)>>6] = _mm512_set1_epi64(q->pw[(i-n)&63]); // Changed from i>>5 to i>>6 and i&31 to i&63
            //print_m512i(q->z_bits[i][i>>6]);
        }
        q->phase[i] = 0;
    }
    std::cout << "Initialized state with " << n << " qubits.\n";
    return;
}

/**
 * Reads a quantum program from a text file and populates the given QProg object.
 * 
 * @param h A pointer to the QProg object to populate.
 * @param fn The filename of the quantum program.
 * @param params Additional parameters for reading the program (optional).
 *              - 'q': Display quantum state after each gate.
 *              - 'p': Display the program after reading.
 *              - 't': Display execution time.
 *              - 's': Silent mode (no output).
 *              - 'm': Suppress measurement results.
 * 
 * @throws std::runtime_error if the file is not found or has an invalid format.
 */
void readprog(QProg *h, std::string fn, std::string params)
{

    
    h->DISPQSTATE = false;
    h->DISPTIME = false;
    h->SILENT = false;
    h->DISPPROG = false;
    h->SUPPRESSM = false;

    if (!params.empty())
    {
        for (int i = 1; i < params.size(); i++)
        {
            if ((tolower(params[i])=='q')) h->DISPQSTATE = true;
            if ((tolower(params[i])=='p')) h->DISPPROG = true;
            if ((tolower(params[i])=='t')) h->DISPTIME = true;
            if ((tolower(params[i])=='s')) h->SILENT = true;
            if ((tolower(params[i])=='m')) h->SUPPRESSM = true;
        }
    }
    
    std::ifstream fp;
    fp.open(fn);
    if (!fp)
    {
        std::ostringstream oss;
        oss << "File not found: " << fn;
        throw std::runtime_error(oss.str());
    }

    std::string line;
    while (std::getline(fp, line) && line != "#");

    if (line != "#") 
    {
        std::ostringstream oss;
        oss << "Invalid file format in file: " << fn;
        throw std::runtime_error(oss.str());
    }

    h->gate_count = 0;
    h->num_qubits = 0;

    while (std::getline(fp, line))
    {
        std::stringstream ss(line);
        char c;
        int val;
        ss >> c;
        if (c == 'c' || c == 'C')
        {
            ss >> val;
            if (val+1 > h->num_qubits) h->num_qubits = val+1;
            ss >> val;
            if (val+1 > h->num_qubits) h->num_qubits = val+1;
        }
        else
        {
            ss >> val;
            if (val+1 > h->num_qubits) h->num_qubits = val+1;
        }
        h->gate_count++;
    }
    fp.close();
    h->op.resize(h->gate_count);
    h->control_qubit.resize(h->gate_count);
    h->target_qubit.resize(h->gate_count);
    fp.open(fn);
    char c;
    while (!fp.eof()&&(c!='#'))
        fp >> c;
    while (!fp.eof())
    {
        std::string gate;
        int qubit1, qubit2;
        fp >> gate;
        if (gate == "c")
        {
            fp >> qubit1 >> qubit2;
            // Store the CNOT gate operation in h
        }
        else if (gate == "h")
        {
            fp >> qubit1;
            // Store the Hadamard gate operation in h
        }
    }
    fp.close();

    return;
}

int main(int argc, char *argv[])
{
    QProg *program;
    QState *qubits;
    int param=0; // whether there are command-line parameters

    std::cout << "\nCHP: Efficient Simulator for Stabilizer Quantum Circuits (multithreaded + AVX-512)";
    std::cout << "\nby Scott Aaronson, modified by Steven Nguyen\n";
    if (argc==1) std::cout<< "\nSyntax: chp [-options] <filename> [input]\n";
    if (argc > 1 && argv[1][0]=='-') param = 1;

    program = new QProg;
    qubits = new QState;
    if (param) readprog(program, argv[2], argv[1]);
    else readprog(program, argv[1],NULL);
    if (argc==(3+param)) initialize_state(qubits, program->num_qubits, argv[2+param]);
    else initialize_state(qubits, program->num_qubits, NULL);

    delete program;
    delete qubits;

    return 0;
}
