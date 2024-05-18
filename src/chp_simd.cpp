#include <omp.h>
#include <chrono> 
#include <iostream> 
#include <string>
#include <vector>
#include <immintrin.h>
#include <bitset>
#include <string.h>
#include <sstream>
#include <fstream>
#include <ctime>
#include "quantum.cpp"
#include "gates.cpp"

#define         CNOT		0
#define         HADAMARD	1
#define         PHASE		2
#define         MEASURE		3


/**
 * Runs the quantum program on the given quantum state.
 * 
 * @param h Pointer to the QProg structure representing the quantum program.
 * @param q Pointer to the QState structure representing the quantum state.
 */
void runprog(struct QProg *h, struct QState *q) {

    int m; // measurement result
    char not_measured = 1;

    auto start_chrono = std::chrono::high_resolution_clock::now(); // Start the chrono timer
    std::time_t start_c = std::time(nullptr); // Start the C timer

    for (int t = 0; t < h->gate_count; t++)
    {
         if (h->op[t]==CNOT) cnot(q,h->control_qubit[t], h->target_qubit[t]);
         if (h->op[t]==HADAMARD) hadamard(q,h->control_qubit[t]);
         if (h->op[t]==PHASE) phase(q,h->control_qubit[t]);
         if (h->op[t]==MEASURE)
         {
                 not_measured = 0;
                 m = measure(q,h->control_qubit[t],h->SUPPRESSM);
                 if (!h->SILENT)
                 {
                    std::cout << "\nOutcome of measuring qubit " << h->control_qubit[t] << ": ";
                    if (m > 1) 
                        std::cout << m - 2 << " (random)";
                    else 
                        std::cout << m;
                 }
         }
    }
    std::cout << std::endl;
    auto end_chrono = std::chrono::high_resolution_clock::now(); // Stop the chrono timer
    std::time_t end_c = std::time(nullptr); // Stop the C timer

    if (h->DISPTIME)
    {
        std::chrono::duration<double> duration_chrono = end_chrono - start_chrono;
        double dt_chrono = duration_chrono.count();
        double dt_c = std::difftime(end_c, start_c);
        std::cout << "\nMeasurement time (chrono): " << dt_chrono << " seconds" << std::endl;
        std::cout << "\nMeasurement time (C): " << dt_c << " seconds" << std::endl;
    }
    return;

}

/**
 * Applies quantum operations to prepare the state of a quantum system based on a given string.
 * The operations applied depend on the characters in the string as follows:
 * - 'Z': Applies a Hadamard gate, followed by two phase gates, and then another Hadamard gate.
 * - 'x': Applies a Hadamard gate.
 * - 'X': Applies a Hadamard gate, followed by two phase gates.
 * - 'y': Applies a Hadamard gate, followed by a phase gate.
 * - 'Y': Applies a Hadamard gate, followed by three phase gates.
 *
 * @param q A pointer to the quantum state structure.
 * @param s The string representing the desired state preparation operations.
 */
void preparestate(struct QState *q, char *s) {
    
	int l = strlen(s);
	for (int b = 0; b < l; b++)
	{
         if (s[b]=='Z')
         {
                 hadamard(q,b);
                 phase(q,b);
                 phase(q,b);
                 hadamard(q,b);
         }
         if (s[b]=='x') hadamard(q,b);
         if (s[b]=='X')
         {
                 hadamard(q,b);
                 phase(q,b);
                 phase(q,b);
         }
         if (s[b]=='y')
         {
                 hadamard(q,b);
                 phase(q,b);
         }
         if (s[b]=='Y')
         {
                 hadamard(q,b);
                 phase(q,b);
                 phase(q,b);
                 phase(q,b);                 
         }
	}
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
    q->over512 = (n>>9) + 1;
    q->x_bits.resize(2 * n + 1);
    q->z_bits.resize(2 * n + 1);
    q->phase.resize(2 * n + 1);
    std::fill(q->phase.begin(), q->phase.end(), 0);
    q->pw[0] = 1;
    for (int i = 1; i < 32; i++)
        q->pw[i] = 2 * q->pw[i - 1];
    for (int i = 0; i < 2 * n + 1; i++) //Iterate through rows of tableau
    {
        unsigned int array[16] = {0, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0};
        q->x_bits[i].resize(q->over512); // Each row needs ceil(n/512) __m512i's to store the qubits
        q->z_bits[i].resize(q->over512);
        for (int j = 0; j < q->over512; j++) // Example: 1024 qubits needs 2 __m512i's, initialize both to all 0's for particular row
        {
            q->x_bits[i][j] = _mm512_setzero_si512();
            q->z_bits[i][j] = _mm512_setzero_si512();
        }
        if (i < n) { // If you are in the destabilizers, set the x bit (top left of tableau)
            int col_index = i >> 9; // Calculate the row index in __m512i
            int int_index = (i & 511) >> 5; // Calculate the bit index within the __m512i
            array[int_index] = q->pw[i & 31];
            q->x_bits[i][col_index] = _mm512_loadu_si512(array);
            array[int_index] = 0;
        }
        else if (i < 2 * n) { // If you are in the stabilizers, set the z bit (bottom right of tableau)
            int col_index = (i - n) >> 9; // Calculate the row index in __m512i
            int int_index = ((i-n) & 511) >>5; // Calculate the bit index within the __m512i
            array[int_index] = q->pw[(i-n) & 31];
            q->z_bits[i][col_index] = _mm512_loadu_si512(array);
            array[int_index] = 0;
        }
        q->phase[i] = 0;
    }
    std::cout << "Initialized state with " << n << " qubits.\n";
    if (s) preparestate(q, s);
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
    std::string directory = "./circuits/";
    std::string fullPath = directory + fn;
    fp.open(fullPath);
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
    int i = 0;
    while (!fp.eof())
    {
        std::string gate;
        int qubit1, qubit2;
        fp >> gate;
        if (gate == "c")
        {
            fp >> qubit1 >> qubit2;
            h->op[i] = CNOT;
            h->control_qubit[i] = qubit1;
            h->target_qubit[i] = qubit2;
        }
        else if (gate == "h")
        {
            fp >> qubit1;
            h->op[i] = HADAMARD;
            h->control_qubit[i] = qubit1;
        }
        else if (gate == "p")
        {
            fp >> qubit1;
            h->op[i] = PHASE;
            h->control_qubit[i] = qubit1;
        }
        else if (gate == "m")
        {
            fp >> qubit1;
            h->op[i] = MEASURE;
            h->control_qubit[i] = qubit1;
        }
        else
        {
            std::ostringstream oss;
            oss << "Invalid gate: " << gate;
            throw std::runtime_error(oss.str());
        }
        i++;
    }
    fp.close();

    return;
}

int main(int argc, char *argv[])
{
    QProg *program;
    QState *qubits;
    int param=0; // whether there are command-line parameters
    srand(time(0));
    std::cout << "\nCHP: Efficient Simulator for Stabilizer Quantum Circuits (multithreaded + AVX-512)";
    std::cout << "\nby Scott Aaronson, modified by Steven Nguyen\n";
    if (argc==1) std::cout<< "\nSyntax: chp [-options] <filename> [input]\n";
    if (argc > 1 && argv[1][0]=='-') param = 1;

    program = new QProg;
    qubits = new QState;
    if (param) 
        readprog(program, argv[2], argv[1]);
    else 
        readprog(program, argv[1],NULL);
    if (argc==(3+param)) 
        initialize_state(qubits, program->num_qubits, argv[2+param]);
    else 
        initialize_state(qubits, program->num_qubits, NULL);
    runprog(program, qubits);

    delete program;
    delete qubits;

    return 0;
}
