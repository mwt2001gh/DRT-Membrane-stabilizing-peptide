//"Towards an RNA/peptides world by the Direct-RNA-Template mechanism：The emergence of membrane-stabilizing peptides in RNA-based protocells"
//by Yu Shi, Chunwu Yu and Wentao Ma*
//C source codes for the simulation program  --- The case is correponding to Fig.4a in the article

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>

/******* RANDOM NUMBER GENERATOR ********/
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
static long ra[M + 1], nd;
void seed(long seed);  // Random number initialization
long randl(long);      // Random number between 0 and parameter
double randd(void);    // Random double between 0 and 1
/****************************************/

#define RECTXT "Cmspg+nsr.txt"    // The recording file
#define LEN sizeof(struct rna)
#define C 2       // Four types of nucleotides
#define G 3
#define A 1
#define U 4
#define LEN_pep sizeof(struct peptide)
#define P 1     // The assumed six types of amino acids in the system
#define Q 2
#define R 3
#define S 4
#define T 5
#define L 6

#define Codon_length 5   // The length of RNA's binding sequence for a specific amino acid
#define PSEQ G,A,C,U,G   // The RNA sequence that could specifically bind the amino acid "P"
#define QSEQ C,G,U,A,C   // The RNA sequence that could specifically bind the amino acid "Q"
#define RSEQ U,G,C,A,U   // The RNA sequence that could specifically bind the amino acid "R"
#define SSEQ U,C,G,A,U   // The RNA sequence that could specifically bind the amino acid "S"
#define TSEQ A,U,G,A,U   // The RNA sequence that could specifically bind the amino acid "T"
#define LSEQ A,G,C,U,A   // The RNA sequence that could specifically bind the amino acid "L"

#define STSEQ SSEQ,TSEQ  // The "RNA gene" sequence that may encode the peptide "ST" (here being assumed to be the membrane-stabilizing peptide, MSP)
#define PRSEQ PSEQ,RSEQ  // The "RNA gene" sequence that may encode the peptide "PR" (here as a control peptide)
#define NSRSEQ C,A,C,C,G,U,A,G,U,A  // The characteristic sequence of a nucleotide synthetase ribozyme (NSR)
#define CTRLSEQ A,C,U,G,A,U,G,U,A,C  // The characteristic sequence of an RNA species without any function (as a control for NSR)

#define STEPNUM 2000000    // Total time steps of the Monte Carlo simulation
#define STAREC 0           // The step to start recording
#define RECINT 1000        // The interval step of recording
#define MAX_RNA_LENGTH 100    // The maximum RNA length allowed in the simulation
#define PEP_LENGTH 2          // The length of the assumed peptides in the simulation
#define MAX_RNA_NUM_IN_GRIDROOM 1000  // The maximum number of RNA allowed in a single grid room

#define INOCUSTEP 10000      // The step for the inoculation of RNA-based protocells
#define INOCU_CELL_NUM 10    // The number of RNA-based protocells which are inoculated
#define INOCU_SEQ_NUM 5      // The number of RNA sequences within the incoculated protocells
#define INOCUSEQ1 STSEQ      // The following four types of RNA sequences are included in the inoculated protocells  
#define INOCUSEQ2 NSRSEQ
#define INOCUSEQ3 PRSEQ
#define INOCUSEQ4 CTRLSEQ

#define SD 5           // The random seed
#define N 30           // The system is defined as an N×N grid 
#define ROOMNUM (N*N)  // The room number in the grid
#define TNPB   50000   // Total nucleotide precursors (quotients in measurement of nucleotides) introduced into the system in the beginning 
#define TAPB   50000   // Total amphiphilic molecule precursors (quotients in measurement of amphiphiles, i.e. membrane components) 
                                                                                             // introduced into the system in the beginning
#define TAAPB  50000   // Total amino acid precursors (quotients in measurement of amino acids) introduced into the system in the beginning 

#define PNF  0.02      // Probability of a nucleotide forming from its precursor (not catalyzed)
#define PNFR 0.5       // Probability of a nucleotide forming from its precursor catalyzed by NSR
#define PAF  0.02      // Probability of an amphiphile forming from its precursor
#define PND  0.05      // Probability of a nucleotide decaying into its precursor
#define PNDE 0.001     // Probability of a nucleotide residue decaying at RNA's chain end
#define PAD  0.01      // Probability of an amphiphile decaying into its precursor 
#define FDW  0.1       // The factor of molecular decay or degradation within the membrane -- for amphiphiles and peptides
#define PRL  0.000001  // Probability of the random ligation of nucleotides and RNAs 
#define PBB  0.00001   // Probability of a phosphodiester bond breaking
#define FDO  20        // The factor of molecular decay or degradation outside protocells -- for nucleotides and RNAs
#define PAT  0.9       // Probability of an RNA template attracting a substrate (nucleotide or oligomer)
#define PFP 0.0001     // Probability of the false base-pairing (relating to RNA's replicating fidelity)
#define PTL 0.5        // Probability of the template-directed ligation of RNA
#define PSP  0.5       // Probability of the separation of a base pair 
#define PMV 0.9        // Probability of the movement of nucleotides,amphiphiles, amino acids and their precursors
#define LAM 200        // The lower limit number of amphiphilic molecules to form a membrane
#define PMF 0.1        // Probability of a membrane forming
#define PALM 0.001     // Probability of an amphiphile leaving the membrane
#define PAJM 0.2       // Probability of an amphiphile joining the membrane
#define PNPP 0.5       // Probability of a nucleotide precursor permeating through the membrane
#define PAPP 0.9       // Probability of an amphiphile precursor permeating through the membrane
#define PCB 0.0002     // Probability of a protocell breaking
#define PCF 0.001      // Probability of two adjacent protocells fusing with each other 
#define PCD 0.05       // Probability of a protocell dividing
#define PMC 0.1        // Probability of a protocell moving
#define RMRW  (pow(p->length1+p->length2+pblength,1/2.0))  // The relationship between the movement of an RNA and its mass

#define PAAF  0.1     // Probability of an amino acid forming from its precursor
#define PAAD  0.2     // Probability of an amino acid decaying into its precursor
#define PAADE 0.1     // Probability of an amino acid residue decaying at a peptide's chain end
#define PPBB  0.01    // Probability of a peptide bond breaking 
#define PAATL 0.5     // Probability of amino acids' ligation on an RNA template (by the DRT mechanism)
#define PPLM  0.1     // Probability of a peptide leaving the membrane 
#define PPJM  0.9     // Probability of a peptide joining the membrane 
#define PAABR 0.9     // Probability of an amino acid binding onto an RNA template
#define PPLR  0.2     // Probability of an amino acid or peptide leaving RNA
#define PAAPP 0.9     // Probability of an amino acid precursor permeating through the membrane 
#define FMSP 1        // The factor concerning the membrane-stabilizing peptide (on preventing amphiphiles from leaving the membrane)
#define PEP_RMRW  (pow(pep->length,1/2.0))  // The relationship between the movement of a peptide and its mass

struct rna           // RNA (a nucleotide can be deemed as an RNA of length 1)
{                   // In the following annotations,"RNA" would include nucleotides, unless being explicitly distinguished.
    char information[2][MAX_RNA_LENGTH];
    int length1;
    int length2;
    int nick[MAX_RNA_LENGTH];
    // When no peptide or amino acid is bound, denoting whether or not adjacent nucleotides on an RNA template are ligated ---  1: no; 0: yes.
    // When a peptide or amino acid is bound，denoting whether or not the location is the binding location ---
    // 2: the start point of the binding (identical to pep->loc); 1: the following points of the binding; 0: no binding.
    struct rna* next;
    struct rna* prior;
    struct peptide* rna_pep_head;   // The linked-list-head of the binding peptides (including amino acids)
};
struct rna* room_head[2][N][N]; // RNA in grid rooms
struct rna* p, * p1, * p2, * p3, * p4, * ps, * ps1, * ps2;
struct rna* rna_arr[MAX_RNA_NUM_IN_GRIDROOM];

struct peptide     // Peptide (an amino acid can be deemed as a peptide of length 1)
{                  // In the following annotations,"peptides" would include amino acids, unless being explicitly distinguished.
    char information[PEP_LENGTH];
    int length;
    int men;   // The status of the peptide --- 2: binding with RNA; 1: within the menbrane;  0: free.
    int loc;   // When men=2, it refers to the start location where it is bound with the RNA.
    struct peptide* next;
    struct peptide* prior;
};
struct peptide* room_pep_head[2][N][N]; // Peptides in grid rooms
struct peptide* pep, * pep1, * pep2, * pep3;

static char p_codseq[50] = { PSEQ };
static char q_codseq[50] = { QSEQ };
static char r_codseq[50] = { RSEQ };
static char s_codseq[50] = { SSEQ };
static char t_codseq[50] = { TSEQ };
static char l_codseq[50] = { LSEQ };
static char pr_codseq[50] = { PRSEQ };
static char st_codseq[50] = { STSEQ };
static char nsr_seq[50] = { NSRSEQ };
static char ctrl_seq[50] = { CTRLSEQ };

static char inocuseq1[50] = { INOCUSEQ1 };
static char inocuseq2[50] = { INOCUSEQ2 };
static char inocuseq3[50] = { INOCUSEQ3 };
static char inocuseq4[50] = { INOCUSEQ4 };

static int pn_arr[2][N][N];  // Precursors of nucleotides in grid rooms
static int pa_arr[2][N][N];  // Precursors of amphiphiles in grid rooms
static int a_arr[2][N][N];   // Amphiphilies in grid rooms
static int m_arr[N][N];      // Amphiphiles within the membrane of protocells in grid rooms
static int paa_arr[2][N][N]; // Precursors of amino acids in grid rooms
int flagcmov[N][N];
int flagcdiv[N][N];

int x, y;                 // The coordinate of rooms in the grid
long step;                // Cycle variable for Monte Carlo steps
int g = 0;				  // Recording times
int over_max_len = 0;
long available;
long availabl[ROOMNUM];
int k, pn_bef, pa_bef, a_bef, m_bef, paa_bef, randcase;
int nsrlength, ctrllength;
int inoculength1, inoculength2, inoculength3, inoculength4;
int p_codlength, q_codlength, r_codlength, s_codlength, t_codlength, l_codlength;
int pq_codlength, pr_codlength, st_codlength, sl_codlength;
int p_loc, q_loc, r_loc, s_loc, t_loc, l_loc;
int flag1, flag2, flag3, flagnsr, flagst, flagpr, flagntsyn, flagamsyn, flagctrl;

// For data recording
float total_nt_mat[(STEPNUM - STAREC) / RECINT + 1];
float total_am_mat[(STEPNUM - STAREC) / RECINT + 1];
float total_pep_num[(STEPNUM - STAREC) / RECINT + 1];
float rna_num[(STEPNUM - STAREC) / RECINT + 1];
float pn_num[(STEPNUM - STAREC) / RECINT + 1];
float am_num[(STEPNUM - STAREC) / RECINT + 1];
float pam_num[(STEPNUM - STAREC) / RECINT + 1];
float paa_num[(STEPNUM - STAREC) / RECINT + 1];
float pep_num[(STEPNUM - STAREC) / RECINT + 1];
float STpepM_num[(STEPNUM - STAREC) / RECINT + 1];
float STpepOM_num[(STEPNUM - STAREC) / RECINT + 1];
float PRpep_num[(STEPNUM - STAREC) / RECINT + 1];
float STpep_rna_num[(STEPNUM - STAREC) / RECINT + 1];
float PRpep_rna_num[(STEPNUM - STAREC) / RECINT + 1];
float nsr_num[(STEPNUM - STAREC) / RECINT + 1];
float ctrl_num[(STEPNUM - STAREC) / RECINT + 1];

float cell_num[(STEPNUM - STAREC) / RECINT + 1];
float cell_nsr_num[(STEPNUM - STAREC) / RECINT + 1];
float cell_st_num[(STEPNUM - STAREC) / RECINT + 1];
float cell_stnsr_num[(STEPNUM - STAREC) / RECINT + 1];
float cell_pr_num[(STEPNUM - STAREC) / RECINT + 1];
float cell_ctrl_num[(STEPNUM - STAREC) / RECINT + 1];

/////////////////////////////////////////////////////////////////////////
void seed(long seed)   // Random number initialization
{
    int a;
    if (seed < 0) { puts("SEED error."); exit(1); }
    ra[0] = (long)fmod(16807.0 * (double)seed, 2147483647.0);
    for (a = 1; a <= M; a++)
    {
        ra[a] = (long)fmod(16807.0 * (double)ra[a - 1], 2147483647.0);
    }
}

/////////////////////////////////////////////////////////////////////////
long randl(long num)      // Random integer number between 0 and num-1
{
    return(RandomInteger % num);
}

/////////////////////////////////////////////////////////////////////////
double randd(void)        // Random real number between 0 and 1
{
    return((double)RandomInteger / RIMAX);
}

/////////////////////////////////////////////////////////////////////////
void avail_xy_init(void)   // Initialization for xy_choose
{
    int j;
    for (j = 0; j < ROOMNUM; j++)
    {
        availabl[j] = j + 1;
    }
    available = ROOMNUM;
}

/////////////////////////////////////////////////////////////////////////
void xy_choose(void)       // Picks a room at random
{
    long rl, s;
    rl = randl(available);
    s = availabl[rl];
    x = (s - 1) % N;
    y = (s - 1) / N;
    availabl[rl] = availabl[available - 1];
    available--;
}

/////////////////////////////////////////////////////////////////////////
void fresh_rna(int h)     // Updating an RNA for the next time step
{
    p1 = p->prior;
    p2 = p->next;
    p3 = room_head[!h][y][x]->next;
    room_head[!h][y][x]->next = p;
    p->next = p3;
    p->prior = room_head[!h][y][x];
    if (p3 != room_head[!h][y][x])p3->prior = p;
    p1->next = p2;
    if (p2 != room_head[h][y][x])p2->prior = p1;
    p = p1;
}

void rna_shuffle(int h) //Random order-arrangement of nodes in the linked list for recording RNA in a grid room
{
    int a, b, len;
    for (a = 0, p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
    {
        if (a == MAX_RNA_NUM_IN_GRIDROOM) { printf("Too many RNA molecules in a single grid room"); exit(0); }
        rna_arr[a] = p;
        a++;
    }
    len = a;
    rna_arr[a] = NULL;
    for (a = len - 1; a >= 1; a--)
    {
        p = rna_arr[a];
        b = RandomInteger % (a);
        rna_arr[a] = rna_arr[b];
        rna_arr[b] = p;
    }
    p = room_head[h][y][x];
    for (a = 0; a < len; a++)
    {
        p->next = rna_arr[a];
        rna_arr[a]->prior = p;
        p = p->next;
    }
    p->next = room_head[h][y][x];
}

void fresh_pep(int h)     // Updating a peptide for the next time step
{
    pep1 = pep->prior;
    pep2 = pep->next;
    pep3 = room_pep_head[!h][y][x]->next;
    room_pep_head[!h][y][x]->next = pep;
    pep->next = pep3;
    pep->prior = room_pep_head[!h][y][x];
    if (pep3 != room_pep_head[!h][y][x])pep3->prior = pep;
    pep1->next = pep2;
    if (pep2 != room_pep_head[h][y][x])pep2->prior = pep1;
    pep = pep1;
}

int findseq_loc(char subseq[], int subseqlength, struct rna* p)  // Find a specific subsequence in an RNA that can bind with an amino acid
{
    int a, b, flag1 = 0, flag2 = 0;

    if (p->length1 >= subseqlength)
    {
        for (b = 0; p->length1 - subseqlength - b >= 0; b++)
        {
            flag2 = 0;
            for (a = 0; a < subseqlength; a++)
            {
                if (p->information[0][b + a] == subseq[a] && p->nick[b + a] == 0)continue; // The locations occupied by amino acids cannot bind them again.
                else { flag2 = 1; break; }
            }
            if (flag2 == 0)break;
        }
        if (flag2 == 1)flag1 = 1;
    }
    else flag1 = 1;

    if (flag1 == 0)return(b);   // Yes, the RNA contains the subsequence,  return the location of the first nucleotide
    else return(-1);   // No, the RNA does not contain the subsequence
}

int findseq_loc2(char subseq[], int subseqlength, struct rna* p)   // Find a specific subsequence in an RNA
{
    int a, b, flag1 = 0, flag2 = 0;

    if (p->length1 >= subseqlength)
    {
        for (b = 0; p->length1 - subseqlength - b >= 0; b++)
        {
            flag2 = 0;
            for (a = 0; a < subseqlength; a++)
            {
                if (p->information[0][b + a] == subseq[a])continue;
                else { flag2 = 1; break; }
            }
            if (flag2 == 0)break;
        }
        if (flag2 == 1)flag1 = 1;
    }
    else flag1 = 1;

    if (flag1 == 0)return(b);   // Yes, the RNA contains the subsequence,  return the location of the first nucleotide
    else return(-1);   // No, the RNA does not contain the subsequence
}

/////////////////////////////////////////////////////////////////////////
void cell_fusing(int ay, int ax)
{
    a_arr[0][(N + y + ay) % N][(N + x + ax) % N] += a_arr[0][y][x]; // Inner amphiphilic molecule moving
    a_arr[0][y][x] = 0;
    pa_arr[0][(N + y + ay) % N][(N + x + ax) % N] += pa_arr[0][y][x]; // Inner amphiphilic molecule precursors moving
    pa_arr[0][y][x] = 0;
    pn_arr[0][(N + y + ay) % N][(N + x + ax) % N] += pn_arr[0][y][x]; // Inner nucleotide precursors moving
    pn_arr[0][y][x] = 0;

    if (room_head[0][y][x]->next != room_head[0][y][x]) // Inner RNA moving
    {
        p = room_head[0][y][x];
        do {
            p = p->next;
        } while (p->next != room_head[0][y][x]); 
        p->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
        if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) p->next->prior = p;
        room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][y][x]->next;
        (room_head[0][y][x]->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_head[0][y][x]->next = room_head[0][y][x];
    }

    paa_arr[0][(N + y + ay) % N][(N + x + ax) % N] += paa_arr[0][y][x]; // Inner amino acid precursors moving
    paa_arr[0][y][x] = 0;

    if (room_pep_head[0][y][x]->next != room_pep_head[0][y][x]) // Inner peptides or peptides within the menbrane moving
    {
        pep = room_pep_head[0][y][x];
        do {
            pep = pep->next;
        } while (pep->next != room_pep_head[0][y][x]); 
        pep->next = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
        if (pep->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]) pep->next->prior = pep;
        room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_pep_head[0][y][x]->next;
        (room_pep_head[0][y][x]->next)->prior = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_pep_head[0][y][x]->next = room_pep_head[0][y][x];
    }

    m_arr[(N + y + ay) % N][(N + x + ax) % N] += m_arr[y][x]; // Membrane fusing
    m_arr[y][x] = 0;
}

/////////////////////////////////////////////////////////////////////////
void oneway_outer_moving(int ay, int by, int ax, int bx)  // Molecules' movement accompanying the cell division or movement in one direction
{
    a_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += a_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer amphiphiles moving
    a_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    pa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += pa_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer amphiphile precusors moving
    pa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    pn_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += pn_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer nucleotide precursors moving
    pn_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

    if (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) // Outer RNA moving
    {
        p = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        do {
            p = p->next;
        } while (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]);
        p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
        if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N])(room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
        room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
        (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
        room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
    }

    paa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N] += paa_arr[0][(N + y + ay) % N][(N + x + ax) % N]; // Outer amino acid precursors moving
    paa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;

    if (room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]) // Outer peptides moving
    {
        pep = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        do {
            pep = pep->next;
        } while (pep->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]);
        pep->next = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
        if (room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N])(room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = pep;
        room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
        (room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
        room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
    }
}

/////////////////////////////////////////////////////////////////////////
void twoway_outer_moving(int ay, int by, int cy, int ax, int bx, int cx) // Molecules' movement accompanying the cell division or movement in two direction
{
    a_bef = a_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amphiphiles moving
    a_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < a_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  a_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  a_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        default: printf("two way outer l moving error in cell division");
        }
    }

    pa_bef = pa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amphiphile precusors moving
    pa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < pa_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  pa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  pa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        default: printf("two way outer pl moving error in cell division");
        }
    }

    pn_bef = pn_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide precursors moving
    pn_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < pn_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  pn_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  pn_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        default: printf("two way outer pn moving error in cell division");
        }
    }

    for (p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; p != room_head[0][(N + y + ay) % N][(N + x + ax) % N]; p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
    {	                                                     // Outer RNA moving
        room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p->next;
        if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (p->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        randcase = randl(2);
        switch (randcase)
        {
        case 0:
            p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
            if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
            room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = p;
            p->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
            break;
        case 1:
            p->next = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
            if (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = p;
            room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = p;
            p->prior = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
            break;
        default: printf("twoway outer moving error");
        }
    }

    paa_bef = paa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amino acid precursors moving
    paa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < paa_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  paa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  paa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        default: printf("two way outer paa moving error in cell division");
        }
    }

    for (pep = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; pep != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]; pep = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
    {	                                                     // Outer peptides moving
        room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = pep->next;
        if (pep->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]) (pep->next)->prior = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        randcase = randl(2);
        switch (randcase)
        {
        case 0:
            pep->next = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
            if (room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = pep;
            room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = pep;
            pep->prior = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
            break;
        case 1:
            pep->next = room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
            if (room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = pep;
            room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = pep;
            pep->prior = room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
            break;
        default: printf("twoway outer moving error");
        }
    }
}

/////////////////////////////////////////////////////////////////////////
void threeway_outer_moving(int ay, int by, int cy, int dy, int ax, int bx, int cx, int dx) // Molecules' movement accompanying the cell division or movement in three direction
{
    a_bef = a_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amphiphiles moving
    a_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < a_bef; k++)
    {
        randcase = randl(3);
        switch (randcase)
        {
        case 0:  a_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  a_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        case 2:  a_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
        default: printf("three way outer l moving error in cell division");
        }
    }

    pa_bef = pa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amphiphile precursors moving
    pa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < pa_bef; k++)
    {
        randcase = randl(3);
        switch (randcase)
        {
        case 0:  pa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  pa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        case 2:  pa_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
        default: printf("three way outer pl moving error in cell division");
        }
    }

    pn_bef = pn_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer nucleotide precursors moving
    pn_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < pn_bef; k++)
    {
        randcase = randl(3);
        switch (randcase)
        {
        case 0:  pn_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  pn_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        case 2:  pn_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
        default: printf("three way outer pn moving error in cell division");
        }
    }

    for (p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; p != room_head[0][(N + y + ay) % N][(N + x + ax) % N]; p = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
    {	                                                    // Outer RNA moving
        room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p->next;
        if (p->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (p->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        randcase = randl(3);
        switch (randcase)
        {
        case 0:
            p->next = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
            if (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = p;
            room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = p;
            p->prior = room_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
            break;
        case 1:
            p->next = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
            if (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = p;
            room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = p;
            p->prior = room_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
            break;
        case 2:
            p->next = room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next;
            if (room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next != room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]) (room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next)->prior = p;
            room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next = p;
            p->prior = room_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N];
            break;
        default: printf("threeway outer moving error");
        }
    }

    paa_bef = paa_arr[0][(N + y + ay) % N][(N + x + ax) % N];   // Outer amino acid precursors moving
    paa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = 0;
    for (k = 0; k < paa_bef; k++)
    {
        randcase = randl(3);
        switch (randcase)
        {
        case 0:  paa_arr[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]++; break;
        case 1:  paa_arr[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]++; break;
        case 2:  paa_arr[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]++; break;
        default: printf("three way outer paa moving error in cell division");
        }
    }

    for (pep = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next; pep != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]; pep = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)
    {	                                                    // Outer peptides moving
        room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = pep->next;
        if (pep->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]) (pep->next)->prior = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        randcase = randl(3);
        switch (randcase)
        {
        case 0:
            pep->next = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next;
            if (room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next != room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]) (room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next)->prior = pep;
            room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N]->next = pep;
            pep->prior = room_pep_head[0][(N + y + ay + by) % N][(N + x + ax + bx) % N];
            break;
        case 1:
            pep->next = room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next;
            if (room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next != room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]) (room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next)->prior = pep;
            room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N]->next = pep;
            pep->prior = room_pep_head[0][(N + y + ay + cy) % N][(N + x + ax + cx) % N];
            break;
        case 2:
            pep->next = room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next;
            if (room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next != room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]) (room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next)->prior = pep;
            room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N]->next = pep;
            pep->prior = room_pep_head[0][(N + y + ay + dy) % N][(N + x + ax + dx) % N];
            break;
        default: printf("threeway outer moving error");
        }
    }
}

/////////////////////////////////////////////////////////////////////////
void cell_moving(int ay, int ax)
{
    a_arr[0][(N + y + ay) % N][(N + x + ax) % N] = a_arr[0][y][x];  // Inner amphiphiles moving
    a_arr[0][y][x] = 0;
    pa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = pa_arr[0][y][x]; // Inner amphiphile precursors moving
    pa_arr[0][y][x] = 0;
    pn_arr[0][(N + y + ay) % N][(N + x + ax) % N] = pn_arr[0][y][x]; // Inner nucleotide precursors moving
    pn_arr[0][y][x] = 0;

    if (room_head[0][y][x]->next != room_head[0][y][x])  // Inner RNA moving
    {
        p = room_head[0][y][x];
        do {
            p = p->next;
        } while (p->next != room_head[0][y][x]); 
        p->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_head[0][y][x]->next;
        (room_head[0][y][x]->next)->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_head[0][y][x]->next = room_head[0][y][x];
    }

    paa_arr[0][(N + y + ay) % N][(N + x + ax) % N] = paa_arr[0][y][x]; // Inner amino acid precursors moving
    paa_arr[0][y][x] = 0;

    if (room_pep_head[0][y][x]->next != room_pep_head[0][y][x])  // Inner peptides or peptides within the menbrane moving
    {
        pep = room_pep_head[0][y][x];
        do {
            pep = pep->next;
        } while (pep->next != room_pep_head[0][y][x]); 
        pep->next = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = room_pep_head[0][y][x]->next;
        (room_pep_head[0][y][x]->next)->prior = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
        room_pep_head[0][y][x]->next = room_pep_head[0][y][x];
    }

    m_arr[(N + y + ay) % N][(N + x + ax) % N] = m_arr[y][x]; // Membrane moving
    m_arr[y][x] = 0;
}

/////////////////////////////////////////////////////////////////////////
void cell_dividing(int ay, int ax)
{
    a_bef = a_arr[0][y][x];  // Inner amphiphiles distributing
    for (k = 0; k < a_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  // Staying in the old cell
            break;
        case 1:  // Moving to the new cell
            a_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; a_arr[0][y][x]--; break;
        default:
            printf("inner am moving error in cell division");
        }
    }

    pa_bef = pa_arr[0][y][x];  // Inner amphiphile precursors distributing
    for (k = 0; k < pa_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  // Staying in the old cell
            break;
        case 1:  // Moving to the new cell
            pa_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; pa_arr[0][y][x]--; break;
        default:
            printf("inner pa moving error in cell division");
        }
    }

    pn_bef = pn_arr[0][y][x];  // Inner nucleotide precursors distributing
    for (k = 0; k < pn_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  // Staying in the old cell
            break;
        case 1:  // Moving to the new cell
            pn_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; pn_arr[0][y][x]--; break;
        default:
            printf("inner pn moving error in cell division");
        }
    }

    for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next) // Inner RNA distributing
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:   // Staying in the old cell
            break;
        case 1:   // Moving to the new cell
            (p->prior)->next = p->next;
            if (p->next != room_head[0][y][x]) (p->next)->prior = p->prior;
            p1 = p;
            p = p->prior;
            p1->next = room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
            if (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_head[0][(N + y + ay) % N][(N + x + ax) % N]) (room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = p1;
            room_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = p1;
            p1->prior = room_head[0][(N + y + ay) % N][(N + x + ax) % N];
            break;
        default: printf("inner RNA moving error");
        }
    }

    paa_bef = paa_arr[0][y][x];  // Inner amino acid precursors distributing
    for (k = 0; k < paa_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  // Staying in the old cell
            break;
        case 1:  // Moving to the new cell
            paa_arr[0][(N + y + ay) % N][(N + x + ax) % N]++; paa_arr[0][y][x]--; break;
        default:
            printf("inner paa moving error in cell division");
        }
    }

    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next) // Inner peptides or peptides within the menbrane distributing
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:   // Staying in the old cell
            break;
        case 1:   // Moving to the new cell
            (pep->prior)->next = pep->next;
            if (pep->next != room_pep_head[0][y][x]) (pep->next)->prior = pep->prior;
            pep1 = pep;
            pep = pep->prior;
            pep1->next = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next;
            if (room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next != room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]) (room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next)->prior = pep1;
            room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N]->next = pep1;
            pep1->prior = room_pep_head[0][(N + y + ay) % N][(N + x + ax) % N];
            break;
        default: printf("inner peptide moving error");
        }
    }

    m_bef = m_arr[y][x] - 2 * LAM;   // Membrane dividing
    m_arr[y][x] = LAM;
    m_arr[(N + y + ay) % N][(N + x + ax) % N] = LAM;
    for (k = 0; k < m_bef; k++)
    {
        randcase = randl(2);
        switch (randcase)
        {
        case 0:  // Staying in the old cell
            m_arr[y][x]++; break;
        case 1:  // Going to the new cell
            m_arr[(N + y + ay) % N][(N + x + ax) % N]++; break;
        default:
            printf("membrane division error in cell division");
        }
    }
}

/////////////////////////////////////////////////////////////////////////
void inits(void)         // Initialization of the system
{
    int j, m, k;
    seed(SD);

    inoculength1 = 0;
    for (j = 0; inocuseq1[j] != 0; j++)
        inoculength1++;
    inoculength2 = 0;
    for (j = 0; inocuseq2[j] != 0; j++)
        inoculength2++;
    inoculength3 = 0;
    for (j = 0; inocuseq3[j] != 0; j++)
        inoculength3++;
    inoculength4 = 0;
    for (j = 0; inocuseq4[j] != 0; j++)
        inoculength4++;

    p_codlength = 0;
    for (j = 0; p_codseq[j] != 0; j++)
        p_codlength++;
    q_codlength = 0;
    for (j = 0; q_codseq[j] != 0; j++)
        q_codlength++;
    r_codlength = 0;
    for (j = 0; r_codseq[j] != 0; j++)
        r_codlength++;
    s_codlength = 0;
    for (j = 0; s_codseq[j] != 0; j++)
        s_codlength++;
    t_codlength = 0;
    for (j = 0; t_codseq[j] != 0; j++)
        t_codlength++;
    l_codlength = 0;
    for (j = 0; l_codseq[j] != 0; j++)
        l_codlength++;
    pr_codlength = 0;
    for (j = 0; pr_codseq[j] != 0; j++)
        pr_codlength++;
    st_codlength = 0;
    for (j = 0; st_codseq[j] != 0; j++)
        st_codlength++;

    for (j = 0; nsr_seq[j] != 0; j++)
        nsrlength++;
    for (j = 0; ctrl_seq[j] != 0; j++)
        ctrllength++;

    for (m = 0; m < 2; m++)     // Initiate the linked-list array of RNA
    {
        for (y = 0; y < N; y++)
        {
            for (x = 0; x < N; x++)
            {
                p1 = (struct rna*)malloc(LEN);
                if (!p1) { printf("\tinit1--memeout\n"); exit(0); }
                room_head[m][y][x] = p1;
                p1->next = room_head[m][y][x];
            }
        }
    }

    for (m = 0; m < 2; m++)     // Initiate the linked-list array of peptides
    {
        for (y = 0; y < N; y++)
        {
            for (x = 0; x < N; x++)
            {
                pep1 = (struct peptide*)malloc(LEN_pep);
                if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                room_pep_head[m][y][x] = pep1;
                pep1->next = room_pep_head[m][y][x];
            }
        }
    }

    for (k = 0; k < TNPB; k++)  // Initial distribution of nucleotide precursors
    {
        x = randl(N);
        y = randl(N);
        pn_arr[0][y][x]++;
    }

    for (k = 0; k < TAPB; k++)  // Initial distribution of amphiphile precursors
    {
        x = randl(N);
        y = randl(N);
        pa_arr[0][y][x]++;
    }

    for (k = 0; k < TAAPB; k++)  // Initial distribution of amino acid precursors
    {
        x = randl(N);
        y = randl(N);
        paa_arr[0][y][x]++;
    }
}

/////////////////////////////////////////////////////////////////////////
void inoculate(void)   // Inoculation of RNA-based protocells
{
    int k0, k1, k2, i;

    for (k0 = 0; k0 < INOCU_CELL_NUM; k0++) {
        x = randl(N);
        y = randl(N);
        for (k1 = 0; k1 < INOCU_SEQ_NUM; k1++) {
            p2 = (struct rna*)malloc(LEN);
            if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
            for (k2 = 0; k2 < inoculength1; k2++) p2->information[0][k2] = inocuseq1[k2];
            p2->information[0][k2] = 0;
            p2->information[1][0] = 0;
            p2->length1 = inoculength1;
            p2->length2 = 0;
            p2->next = room_head[0][y][x]->next;
            if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
            room_head[0][y][x]->next = p2;
            p2->prior = room_head[0][y][x];

            for (i = 0; i < MAX_RNA_LENGTH; i++) p2->nick[i] = 0;

            pep1 = (struct peptide*)malloc(LEN_pep);
            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
            p2->rna_pep_head = pep1;
            pep1->next = p2->rna_pep_head;

            p2 = (struct rna*)malloc(LEN);
            if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
            for (k2 = 0; k2 < inoculength2; k2++) p2->information[0][k2] = inocuseq2[k2];
            p2->information[0][k2] = 0;
            p2->information[1][0] = 0;
            p2->length1 = inoculength2;
            p2->length2 = 0;
            p2->next = room_head[0][y][x]->next;
            if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
            room_head[0][y][x]->next = p2;
            p2->prior = room_head[0][y][x];

            for (i = 0; i < MAX_RNA_LENGTH; i++) p2->nick[i] = 0;

            pep1 = (struct peptide*)malloc(LEN_pep);
            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
            p2->rna_pep_head = pep1;
            pep1->next = p2->rna_pep_head;

            p2 = (struct rna*)malloc(LEN);
            if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
            for (k2 = 0; k2 < inoculength3; k2++) p2->information[0][k2] = inocuseq3[k2];
            p2->information[0][k2] = 0;
            p2->information[1][0] = 0;
            p2->length1 = inoculength3;
            p2->length2 = 0;
            p2->next = room_head[0][y][x]->next;
            if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
            room_head[0][y][x]->next = p2;
            p2->prior = room_head[0][y][x];

            for (i = 0; i < MAX_RNA_LENGTH; i++) p2->nick[i] = 0;

            pep1 = (struct peptide*)malloc(LEN_pep);
            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
            p2->rna_pep_head = pep1;
            pep1->next = p2->rna_pep_head;

            p2 = (struct rna*)malloc(LEN);
            if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
            for (k2 = 0; k2 < inoculength4; k2++) p2->information[0][k2] = inocuseq4[k2];
            p2->information[0][k2] = 0;
            p2->information[1][0] = 0;
            p2->length1 = inoculength4;
            p2->length2 = 0;
            p2->next = room_head[0][y][x]->next;
            if (p2->next != room_head[0][y][x])(p2->next)->prior = p2;
            room_head[0][y][x]->next = p2;
            p2->prior = room_head[0][y][x];

            for (i = 0; i < MAX_RNA_LENGTH; i++) p2->nick[i] = 0;

            pep1 = (struct peptide*)malloc(LEN_pep);
            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
            p2->rna_pep_head = pep1;
            pep1->next = p2->rna_pep_head;
        }
        m_arr[y][x] += LAM * 1.5;
    }
    printf("\ninoculating\n");
}

////////////////////////////////////////////////////////////////////////////////////////
void unit_action(void)      // Action (movement and events) of units (molecules and protocells) in the system
{
    int a, b, d, i, j, randnt, nt_turn, left, right, up, down, tmr_xy, tmpep_xy, pblength;
    double f, f1, rtdaddlig, rtdaddphili, rntsyn, ramsyn;

    //=============================== Movement of molecules
    //-----------------------------------Amphiphiles moving (including joining into the membrane)
    for (y = 0; y < N; y++)   // Considering each room in the grid
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
            {					//  Possibly joining into the membrane
                a_bef = a_arr[0][y][x];
                a_arr[0][y][x] = 0;
                for (k = 0; k < a_bef; k++)
                {
                    if (randd() < PMV && randd() < PAJM) m_arr[y][x]++;
                    else a_arr[1][y][x]++;
                }
            }
            else   // No membrane
            {		// Possibly moving, possibly joining into membrane of a protocell at an adjacent room
                a_bef = a_arr[0][y][x];
                a_arr[0][y][x] = 0;
                for (k = 0; k < a_bef; k++)
                {
                    if (randd() < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:   // To left
                            if (m_arr[y][(N + x - 1) % N] == 0) a_arr[1][y][(N + x - 1) % N]++; // Moving
                            else
                            {
                                if (randd() < PAJM) m_arr[y][(N + x - 1) % N]++; // Joining
                                else a_arr[1][y][x]++;
                            }
                            break;
                        case 1:  // To right
                            if (m_arr[y][(x + 1) % N] == 0) a_arr[1][y][(x + 1) % N]++;
                            else
                            {
                                if (randd() < PAJM) m_arr[y][(x + 1) % N]++;
                                else a_arr[1][y][x]++;
                            }
                            break;
                        case 2:  // To up
                            if (m_arr[(N + y - 1) % N][x] == 0) a_arr[1][(N + y - 1) % N][x]++;
                            else
                            {
                                if (randd() < PAJM) m_arr[(N + y - 1) % N][x]++;
                                else a_arr[1][y][x]++;
                            }
                            break;
                        case 3:   // To down
                            if (m_arr[(y + 1) % N][x] == 0) a_arr[1][(y + 1) % N][x]++;
                            else
                            {
                                if (randd() < PAJM) m_arr[(y + 1) % N][x]++;
                                else a_arr[1][y][x]++;
                            }
                            break;
                        default: printf("l moving error");
                        }
                    }
                    else a_arr[1][y][x]++;
                }
            }
        }
    }

    //-------------------------------- Amphiphiles leaving the membrane
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] == 0)continue;  // No membrane
            tmr_xy = 0;
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
                tmr_xy += p->length1 + p->length2; // Caculating total materials of RNA in the protocell
            // To consider the effect of osomitic pressure caused by the impermeable materials in the protocell
            tmpep_xy = 0;
            for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
            {
                if (pep->men == 1 && pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))
                    tmpep_xy++;   // Caculating the number of MSP within the membrane
            }

            m_bef = m_arr[y][x];
            for (k = 0; k < m_bef; k++)
            {
                if (randd() < (PALM / (1 + tmr_xy / pow(m_arr[y][x] / 2.0, 1.5))) / (1 + FMSP * tmpep_xy))
                {     // The more MSP within the membrane, the more difficult for amphiphiles to leave the membrane.
                    m_arr[y][x]--;  // Leaving
                    if (randd() < 0.5) a_arr[1][y][x]++;  // Going into the protocell
                    else                              // Going out of the protocell
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) a_arr[1][y][(N + x - 1) % N]++;
                            else m_arr[y][x]++;
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) a_arr[1][y][(x + 1) % N]++;
                            else m_arr[y][x]++;
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) a_arr[1][(N + y - 1) % N][x]++;
                            else m_arr[y][x]++;
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) a_arr[1][(y + 1) % N][x]++;
                            else m_arr[y][x]++;
                            break;
                        default: printf("amphiphilic molecule leaving error");
                        }
                    }
                }
            }
        }
    }

    //---------------------------------- Peptides moving (including joining into the membrane)
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0)  // Membrane existing (the room is occupied by a protocell)
            {
                for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                {
                    if (pep->men != 0)  // Only free peptides are considered below
                    {
                        fresh_pep(0);
                        continue;
                    }

                    if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))  // Only MSP can join the membrane.
                    {
                        if (randd() < PPJM)
                        {
                            pep->men = 1;
                            fresh_pep(0);
                        }
                        else fresh_pep(0);

                    }
                    else fresh_pep(0);

                }
            }
            else   // No membrane
            {		// Possibly moving
                for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                {
                    if (pep->men != 0)  // Only free peptides are considered below.
                    {
                        fresh_pep(0);
                        continue;
                    }
                    if (randd() * PEP_RMRW < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[1][y][(N + x - 1) % N]->next;
                                room_pep_head[1][y][(N + x - 1) % N]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[1][y][(N + x - 1) % N];
                                if (pep3 != room_pep_head[1][y][(N + x - 1) % N])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))  // Only the membrane-stabilizing peptide can enter the membrane
                            {
                                if (randd() < PPJM)
                                {
                                    pep->men = 1;
                                    pep1 = pep->prior;
                                    pep2 = pep->next;
                                    pep3 = room_pep_head[1][y][(N + x - 1) % N]->next;
                                    room_pep_head[1][y][(N + x - 1) % N]->next = pep;
                                    pep->next = pep3;
                                    pep->prior = room_pep_head[1][y][(N + x - 1) % N];
                                    if (pep3 != room_pep_head[1][y][(N + x - 1) % N])pep3->prior = pep;
                                    pep1->next = pep2;
                                    if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                    pep = pep1;
                                }
                                else fresh_pep(0);
                            }
                            else fresh_pep(0);
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[1][y][(x + 1) % N]->next;
                                room_pep_head[1][y][(x + 1) % N]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[1][y][(x + 1) % N];
                                if (pep3 != room_pep_head[1][y][(x + 1) % N])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))
                            {
                                if (randd() < PPJM)
                                {
                                    pep->men = 1;

                                    pep1 = pep->prior;
                                    pep2 = pep->next;
                                    pep3 = room_pep_head[1][y][(x + 1) % N]->next;
                                    room_pep_head[1][y][(x + 1) % N]->next = pep;
                                    pep->next = pep3;
                                    pep->prior = room_pep_head[1][y][(x + 1) % N];
                                    if (pep3 != room_pep_head[1][y][(x + 1) % N])pep3->prior = pep;
                                    pep1->next = pep2;
                                    if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                    pep = pep1;
                                }
                                else fresh_pep(0);

                            }
                            else fresh_pep(0);
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[1][(N + y - 1) % N][x]->next;
                                room_pep_head[1][(N + y - 1) % N][x]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[1][(N + y - 1) % N][x];
                                if (pep3 != room_pep_head[1][(N + y - 1) % N][x])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))
                            {
                                if (randd() < PPJM)
                                {
                                    pep->men = 1;

                                    pep1 = pep->prior;
                                    pep2 = pep->next;
                                    pep3 = room_pep_head[1][(N + y - 1) % N][x]->next;
                                    room_pep_head[1][(N + y - 1) % N][x]->next = pep;
                                    pep->next = pep3;
                                    pep->prior = room_pep_head[1][(N + y - 1) % N][x];
                                    if (pep3 != room_pep_head[1][(N + y - 1) % N][x])pep3->prior = pep;
                                    pep1->next = pep2;
                                    if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                    pep = pep1;

                                }
                                else fresh_pep(0);

                            }
                            else fresh_pep(0);
                            break;

                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[1][(y + 1) % N][x]->next;
                                room_pep_head[1][(y + 1) % N][x]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[1][(y + 1) % N][x];
                                if (pep3 != room_pep_head[1][(y + 1) % N][x])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))
                            {
                                if (randd() < PPJM)
                                {
                                    pep->men = 1;

                                    pep1 = pep->prior;
                                    pep2 = pep->next;
                                    pep3 = room_pep_head[1][(y + 1) % N][x]->next;
                                    room_pep_head[1][(y + 1) % N][x]->next = pep;
                                    pep->next = pep3;
                                    pep->prior = room_pep_head[1][(y + 1) % N][x];
                                    if (pep3 != room_pep_head[1][(y + 1) % N][x])pep3->prior = pep;
                                    pep1->next = pep2;
                                    if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;
                                    pep = pep1;
                                }
                                else fresh_pep(0);
                            }
                            else fresh_pep(0);
                            break;

                        default: printf("rna moving error");
                        }
                    }
                    else fresh_pep(0);
                }
            }
        }
    }

    // ------------------------------ Peptides leaving the membrane
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] == 0)  // No membrane
            {
                for (pep = room_pep_head[1][y][x]->next; pep != room_pep_head[1][y][x]; pep = pep->next)
                {
                    fresh_pep(1);
                }
                continue;
            }
            tmr_xy = 0;
            for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
                tmr_xy += p->length1 + p->length2; // Caculating total materials of RNA in the protocell
            // To consider the effect of osomitic pressure caused by the impermeable materials in the protocell

            for (pep = room_pep_head[1][y][x]->next; pep != room_pep_head[1][y][x]; pep = pep->next)
            {
                if (pep->men != 1)  // Only the peptides within the membrane are considered below.
                {
                    fresh_pep(1);
                    continue;
                }

                if (randd() < PPLM / (1 + tmr_xy / pow(m_arr[y][x] / 2.0, 1.5)))
                {
                    pep->men = 0;  // Leaving
                    if (randd() < 0.5) fresh_pep(1);  // Going into the protocell
                    else                              // Going out of the protocell
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[0][y][(N + x - 1) % N]->next;
                                room_pep_head[0][y][(N + x - 1) % N]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[0][y][(N + x - 1) % N];
                                if (pep3 != room_pep_head[0][y][(N + x - 1) % N])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[1][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else   // Staying within the membrane
                            {
                                pep->men = 1;
                                fresh_pep(1);
                            }
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[0][y][(x + 1) % N]->next;
                                room_pep_head[0][y][(x + 1) % N]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[0][y][(x + 1) % N];
                                if (pep3 != room_pep_head[0][y][(x + 1) % N])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[1][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else
                            {
                                pep->men = 1;
                                fresh_pep(1);
                            }
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[0][(N + y - 1) % N][x]->next;
                                room_pep_head[0][(N + y - 1) % N][x]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[0][(N + y - 1) % N][x];
                                if (pep3 != room_pep_head[0][(N + y - 1) % N][x])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[1][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else
                            {
                                pep->men = 1;
                                fresh_pep(1);
                            }
                            break;

                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0)
                            {
                                pep1 = pep->prior;
                                pep2 = pep->next;
                                pep3 = room_pep_head[0][(y + 1) % N][x]->next;
                                room_pep_head[0][(y + 1) % N][x]->next = pep;
                                pep->next = pep3;
                                pep->prior = room_pep_head[0][(y + 1) % N][x];
                                if (pep3 != room_pep_head[0][(y + 1) % N][x])pep3->prior = pep;
                                pep1->next = pep2;
                                if (pep2 != room_pep_head[1][y][x])pep2->prior = pep1;
                                pep = pep1;
                            }
                            else
                            {
                                pep->men = 1;
                                fresh_pep(1);
                            }
                            break;

                        default: printf("peptides leaving error");
                        }
                    }
                }
                else fresh_pep(1);
            }
        }
    }

    //----------------------------- Binding of amino acids with RNA
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
            {

                if (p->length2 != 0) continue;  // Only single-chain RNAs are considered below

                q_loc = findseq_loc(q_codseq, q_codlength, p);  // Considering the binding of amino acid Q onto RNA
                if (q_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue; // Only amino acids are considered below 
                        if (pep->men != 0) continue;   // Only free amino acids are considered below
                        if (pep->information[0] != 2) continue;  // Only Q is considered below
                        if (randd() < PAABR)
                        {
                            p->nick[q_loc] = 2;
                            for (i = q_loc + 1; i < q_loc + q_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;  // Binding with RNA
                            pep->loc = q_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)  // RNA without peptides
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else    // For RNA already binding with peptides
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > q_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }

                r_loc = findseq_loc(r_codseq, r_codlength, p); // Considering the binding of amino acid R onto RNA
                if (r_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue;
                        if (pep->men != 0) continue;
                        if (pep->information[0] != 3) continue;
                        if (randd() < PAABR)
                        {
                            p->nick[r_loc] = 2;
                            for (i = r_loc + 1; i < r_loc + r_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;
                            pep->loc = r_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > r_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }

                p_loc = findseq_loc(p_codseq, p_codlength, p); // Considering the binding of amino acid P onto RNA
                if (p_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue;
                        if (pep->men != 0) continue;
                        if (pep->information[0] != 1) continue;
                        if (randd() < PAABR)
                        {
                            p->nick[p_loc] = 2;
                            for (i = p_loc + 1; i < p_loc + p_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;
                            pep->loc = p_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            //
                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > p_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }

                s_loc = findseq_loc(s_codseq, s_codlength, p);  // Considering the binding of amino acid S onto RNA
                if (s_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue;
                        if (pep->men != 0) continue;
                        if (pep->information[0] != 4) continue;
                        if (randd() < PAABR)
                        {
                            p->nick[s_loc] = 2;
                            for (i = s_loc + 1; i < s_loc + s_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;
                            pep->loc = s_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > s_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }

                t_loc = findseq_loc(t_codseq, t_codlength, p);   // Considering the binding of amino acid T onto RNA
                if (t_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue;
                        if (pep->men != 0) continue;
                        if (pep->information[0] != 5) continue;
                        if (randd() < PAABR)
                        {
                            p->nick[t_loc] = 2;
                            for (i = t_loc + 1; i < t_loc + t_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;
                            pep->loc = t_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > t_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }

                l_loc = findseq_loc(l_codseq, l_codlength, p);   // Considering the binding of amino acid L onto RNA
                if (l_loc != -1)
                {
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                    {
                        if (pep->length != 1) continue;
                        if (pep->men != 0) continue;
                        if (pep->information[0] != 6) continue;
                        if (randd() < PAABR)
                        {
                            p->nick[l_loc] = 2;
                            for (i = l_loc + 1; i < l_loc + l_codlength; i++) p->nick[i] = 1;
                            pep->men = 2;
                            pep->loc = l_loc;

                            pep1 = pep->prior;
                            pep2 = pep->next;
                            pep1->next = pep2;
                            if (pep2 != room_pep_head[0][y][x])pep2->prior = pep1;

                            pep3 = pep;
                            pep = pep->prior;

                            if (p->rna_pep_head->next == p->rna_pep_head)
                            {
                                p->rna_pep_head->next = pep3;
                                pep3->prior = p->rna_pep_head;
                                pep3->next = p->rna_pep_head;
                            }
                            else
                            {
                                for (pep1 = p->rna_pep_head->next; pep1 != p->rna_pep_head; pep1 = pep1->next)
                                {
                                    if (pep1->loc > l_loc) break;
                                }
                                if (pep1 == p->rna_pep_head)
                                {
                                    pep1 = p->rna_pep_head->next;
                                    while (pep1->next != p->rna_pep_head) pep1 = pep1->next;
                                    pep3->next = pep1->next;
                                    pep1->next = pep3;
                                    pep3->prior = pep1;
                                }
                                else
                                {
                                    pep1->prior->next = pep3;
                                    pep3->prior = pep1->prior;
                                    pep1->prior = pep3;
                                    pep3->next = pep1;
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    //---------------------------------- Adjacent amino acids on the RNA template ligating to form peptides
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
            {
                if (p->length2 != 0) continue;  // Only considering single RNA chains below
                if (p->rna_pep_head->next == p->rna_pep_head) continue; // Only considering RNA bound with amino acids or peptides
                for (pep = p->rna_pep_head->next; pep != p->rna_pep_head; pep = pep->next)
                {
                    if (pep->next == p->rna_pep_head) continue;  // At the end of the peptide chain
                    if (pep->length != 1 || pep->next->length != 1) continue;   // Only single amino acids can be ligated.
                    if (pep->next->loc - pep->loc != Codon_length) continue;   // Only adjacent amino acids can be ligated.
                    if (randd() < PAATL)
                    {
                        pep1 = pep->next;
                        p->nick[pep1->loc] = 1;
                        pep->information[1] = pep1->information[0];
                        pep->length++;
                        pep->next = pep1->next;
                        if (pep1->next != p->rna_pep_head) pep1->next->prior = pep;
                        free(pep1);
                    }
                }
            }
        }
    }

    //-------------------- The separation of peptides and RNA 
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
            {
                if (p->length2 != 0) continue;
                if ((p->rna_pep_head)->next == p->rna_pep_head) continue;
                for (pep = p->rna_pep_head->next; pep != p->rna_pep_head; pep = pep->next)
                {
                    if (randd() < PPLR)  // No "length" factor is added considering the formation of peptides would weaken the amino acid-RNA binding
                    {
                        for (k = pep->loc; k < pep->loc + pep->length * Codon_length; k++) p->nick[k] = 0;
                        pep->men = 0;
                        pep1 = pep->prior;
                        pep2 = pep->next;
                        pep3 = room_pep_head[0][y][x]->next;
                        room_pep_head[0][y][x]->next = pep;
                        pep->next = pep3;
                        pep->prior = room_pep_head[0][y][x];
                        if (pep3 != room_pep_head[0][y][x])pep3->prior = pep;
                        pep1->next = pep2;
                        if (pep2 != p->rna_pep_head)pep2->prior = pep1;
                        pep = pep1;
                    }
                }
            }
        }
    }

    //------------------------------ RNA moving
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0) // Membrane existing
            {					// Not moving
                if (room_head[0][y][x]->next != room_head[0][y][x])
                {
                    p = room_head[0][y][x];
                    do {
                        p = p->next;
                    } while (p->next != room_head[0][y][x]);
                    p->next = room_head[1][y][x];
                    room_head[1][y][x]->next = room_head[0][y][x]->next;
                    (room_head[0][y][x]->next)->prior = room_head[1][y][x];
                    room_head[0][y][x]->next = room_head[0][y][x];
                }
            }
            else   // No membrane
            {		// Possibly moving
                for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
                {
                    pblength = 0;
                    for (pep = p->rna_pep_head->next; pep != p->rna_pep_head; pep = pep->next) //Caculating the mass of amino acids and peptides bound on the RNA
                    {
                        pblength += pep->length;
                    }
                    if (randd() * RMRW < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0)
                            {
                                p1 = p->prior;
                                p2 = p->next;
                                p3 = room_head[1][y][(N + x - 1) % N]->next;
                                room_head[1][y][(N + x - 1) % N]->next = p;
                                p->next = p3;
                                p->prior = room_head[1][y][(N + x - 1) % N];
                                if (p3 != room_head[1][y][(N + x - 1) % N])p3->prior = p;
                                p1->next = p2;
                                if (p2 != room_head[0][y][x])p2->prior = p1;
                                p = p1;
                            }
                            else fresh_rna(0);
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0)
                            {
                                p1 = p->prior;
                                p2 = p->next;
                                p3 = room_head[1][y][(x + 1) % N]->next;
                                room_head[1][y][(x + 1) % N]->next = p;
                                p->next = p3;
                                p->prior = room_head[1][y][(x + 1) % N];
                                if (p3 != room_head[1][y][(x + 1) % N])p3->prior = p;
                                p1->next = p2;
                                if (p2 != room_head[0][y][x])p2->prior = p1;
                                p = p1;
                            }
                            else fresh_rna(0);
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0)
                            {
                                p1 = p->prior;
                                p2 = p->next;
                                p3 = room_head[1][(N + y - 1) % N][x]->next;
                                room_head[1][(N + y - 1) % N][x]->next = p;
                                p->next = p3;
                                p->prior = room_head[1][(N + y - 1) % N][x];
                                if (p3 != room_head[1][(N + y - 1) % N][x])p3->prior = p;
                                p1->next = p2;
                                if (p2 != room_head[0][y][x])p2->prior = p1;
                                p = p1;
                            }
                            else fresh_rna(0);
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0)
                            {
                                p1 = p->prior;
                                p2 = p->next;
                                p3 = room_head[1][(y + 1) % N][x]->next;
                                room_head[1][(y + 1) % N][x]->next = p;
                                p->next = p3;
                                p->prior = room_head[1][(y + 1) % N][x];
                                if (p3 != room_head[1][(y + 1) % N][x])p3->prior = p;
                                p1->next = p2;
                                if (p2 != room_head[0][y][x])p2->prior = p1;
                                p = p1;
                            }
                            else fresh_rna(0);
                            break;
                        default: printf("rna moving error");
                        }
                    }
                    else fresh_rna(0);
                }
            }
        }
    }

    //-------------------------- Amphiphilic molecule precursors moving
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0)  // Membrane existing
            {					// Possibly moving and permeating out
                pa_bef = pa_arr[0][y][x];
                pa_arr[0][y][x] = 0;
                for (k = 0; k < pa_bef; k++)
                {
                    if (randd() < PMV && randd() < PAPP)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) pa_arr[1][y][(N + x - 1) % N]++;
                            else pa_arr[1][y][x]++;
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) pa_arr[1][y][(x + 1) % N]++;
                            else pa_arr[1][y][x]++;
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) pa_arr[1][(N + y - 1) % N][x]++;
                            else pa_arr[1][y][x]++;
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) pa_arr[1][(y + 1) % N][x]++;
                            else pa_arr[1][y][x]++;
                            break;
                        default: printf("pa moving error");
                        }
                    }
                    else pa_arr[1][y][x]++;
                }
            }
            else    // no membrane
            {		// Possibly moving, possibly permeating into a protocell at an adjacent room
                pa_bef = pa_arr[0][y][x];
                pa_arr[0][y][x] = 0;
                for (k = 0; k < pa_bef; k++)
                {
                    if (randd() < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) pa_arr[1][y][(N + x - 1) % N]++;
                            else
                            {
                                if (randd() < PAPP * (m_arr[y][(N + x - 1) % N] / LAM))
                                    pa_arr[1][y][(N + x - 1) % N]++;
                                else pa_arr[1][y][x]++;
                            }
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) pa_arr[1][y][(x + 1) % N]++;
                            else
                            {
                                if (randd() < PAPP * (m_arr[y][(x + 1) % N] / LAM))
                                    pa_arr[1][y][(x + 1) % N]++;
                                else pa_arr[1][y][x]++;
                            }
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) pa_arr[1][(N + y - 1) % N][x]++;
                            else
                            {
                                if (randd() < PAPP * (m_arr[(N + y - 1) % N][x] / LAM))
                                    pa_arr[1][(N + y - 1) % N][x]++;
                                else pa_arr[1][y][x]++;
                            }
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) pa_arr[1][(y + 1) % N][x]++;
                            else
                            {
                                if (randd() < PAPP * (m_arr[(y + 1) % N][x] / LAM))
                                    pa_arr[1][(y + 1) % N][x]++;
                                else pa_arr[1][y][x]++;
                            }
                            break;
                        default: printf("pl moving error");
                        }
                    }
                    else pa_arr[1][y][x]++;
                }
            }
        }
    }

    //--------------------------------- Nucleotide precursors moving
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0)  // Membrane existing
            {				   // Possibly moving and permeating out
                pn_bef = pn_arr[0][y][x];
                pn_arr[0][y][x] = 0;
                for (k = 0; k < pn_bef; k++)
                {
                    if (randd() < PMV && randd() < PNPP)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) pn_arr[1][y][(N + x - 1) % N]++;
                            else pn_arr[1][y][x]++;
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) pn_arr[1][y][(x + 1) % N]++;
                            else pn_arr[1][y][x]++;
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) pn_arr[1][(N + y - 1) % N][x]++;
                            else pn_arr[1][y][x]++;
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) pn_arr[1][(y + 1) % N][x]++;
                            else pn_arr[1][y][x]++;
                            break;
                        default: printf("pn moving error");
                        }
                    }
                    else pn_arr[1][y][x]++;
                }
            }
            else    // no membrane
            {		// Possible moving, possibly permeating into a protocell at an adjacent room
                pn_bef = pn_arr[0][y][x];
                pn_arr[0][y][x] = 0;
                for (k = 0; k < pn_bef; k++)
                {
                    if (randd() < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) pn_arr[1][y][(N + x - 1) % N]++; // Moving
                            else
                            {
                                tmr_xy = 0;
                                for (p = room_head[1][y][(N + x - 1) % N]->next; p != room_head[1][y][(N + x - 1) % N]; p = p->next)
                                    tmr_xy += p->length1 + p->length2; // Caculating total materials of nucleotides and RNA in the target protocell
                                // To consider the effect of Donnan's equilibrium
                                if (randd() < PNPP * (m_arr[y][(N + x - 1) % N] / LAM) / (1.0 + tmr_xy / pow(m_arr[y][(N + x - 1) % N] / 2.0, 1.5)))
                                    pn_arr[1][y][(N + x - 1) % N]++;
                                else pn_arr[1][y][x]++;
                            }
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) pn_arr[1][y][(x + 1) % N]++;
                            else
                            {
                                tmr_xy = 0;
                                for (p = room_head[1][y][(x + 1) % N]->next; p != room_head[1][y][(x + 1) % N]; p = p->next)
                                    tmr_xy += p->length1 + p->length2;
                                if (randd() < PNPP * (m_arr[y][(x + 1) % N] / LAM) / (1.0 + tmr_xy / pow(m_arr[y][(x + 1) % N] / 2.0, 1.5)))
                                    pn_arr[1][y][(x + 1) % N]++;
                                else pn_arr[1][y][x]++;
                            }
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) pn_arr[1][(N + y - 1) % N][x]++;
                            else
                            {
                                tmr_xy = 0;
                                for (p = room_head[1][(N + y - 1) % N][x]->next; p != room_head[1][(N + y - 1) % N][x]; p = p->next)
                                    tmr_xy += p->length1 + p->length2;
                                if (randd() < PNPP * (m_arr[(N + y - 1) % N][x] / LAM) / (1.0 + tmr_xy / pow(m_arr[(N + y - 1) % N][x] / 2.0, 1.5)))
                                    pn_arr[1][(N + y - 1) % N][x]++;
                                else pn_arr[1][y][x]++;
                            }
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) pn_arr[1][(y + 1) % N][x]++;
                            else
                            {
                                tmr_xy = 0;
                                for (p = room_head[1][(y + 1) % N][x]->next; p != room_head[1][(y + 1) % N][x]; p = p->next)
                                    tmr_xy += p->length1 + p->length2;
                                if (randd() < PNPP * (m_arr[(y + 1) % N][x] / LAM) / (1.0 + tmr_xy / pow(m_arr[(y + 1) % N][x] / 2.0, 1.5)))
                                    pn_arr[1][(y + 1) % N][x]++;
                                else pn_arr[1][y][x]++;
                            }
                            break;
                        default: printf("pn moving error");
                        }
                    }
                    else pn_arr[1][y][x]++;
                }
            }
        }
    }

    //-------------------------------- Amino acid precursors moving
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] > 0)  // Membrane existing
            {					// Possibly moving and permeating out
                paa_bef = paa_arr[0][y][x];
                paa_arr[0][y][x] = 0;
                for (k = 0; k < paa_bef; k++)
                {
                    if (randd() < PMV && randd() < PAAPP)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) paa_arr[1][y][(N + x - 1) % N]++;
                            else paa_arr[1][y][x]++;
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) paa_arr[1][y][(x + 1) % N]++;
                            else paa_arr[1][y][x]++;
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) paa_arr[1][(N + y - 1) % N][x]++;
                            else paa_arr[1][y][x]++;
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) paa_arr[1][(y + 1) % N][x]++;
                            else paa_arr[1][y][x]++;
                            break;
                        default: printf("pa moving error");
                        }
                    }
                    else paa_arr[1][y][x]++;
                }
            }
            else    // no membrane
            {		// Possibly moving, possibly permeating into a protocell at an adjacent room
                paa_bef = paa_arr[0][y][x];
                paa_arr[0][y][x] = 0;
                for (k = 0; k < paa_bef; k++)
                {
                    if (randd() < PMV)
                    {
                        randcase = randl(4);   // Four possible directions
                        switch (randcase)
                        {
                        case 0:
                            if (m_arr[y][(N + x - 1) % N] == 0) paa_arr[1][y][(N + x - 1) % N]++;
                            else
                            {
                                if (randd() < PAAPP * (m_arr[y][(N + x - 1) % N] / LAM))
                                    paa_arr[1][y][(N + x - 1) % N]++;
                                else paa_arr[1][y][x]++;
                            }
                            break;
                        case 1:
                            if (m_arr[y][(x + 1) % N] == 0) paa_arr[1][y][(x + 1) % N]++;
                            else
                            {
                                if (randd() < PAAPP * (m_arr[y][(x + 1) % N] / LAM))
                                    paa_arr[1][y][(x + 1) % N]++;
                                else paa_arr[1][y][x]++;
                            }
                            break;
                        case 2:
                            if (m_arr[(N + y - 1) % N][x] == 0) paa_arr[1][(N + y - 1) % N][x]++;
                            else
                            {
                                if (randd() < PAAPP * (m_arr[(N + y - 1) % N][x] / LAM))
                                    paa_arr[1][(N + y - 1) % N][x]++;
                                else paa_arr[1][y][x]++;
                            }
                            break;
                        case 3:
                            if (m_arr[(y + 1) % N][x] == 0) paa_arr[1][(y + 1) % N][x]++;
                            else
                            {
                                if (randd() < PAAPP * (m_arr[(y + 1) % N][x] / LAM))
                                    paa_arr[1][(y + 1) % N][x]++;
                                else paa_arr[1][y][x]++;
                            }
                            break;
                        default: printf("pl moving error");
                        }
                    }
                    else paa_arr[1][y][x]++;
                }
            }
        }
    }
    //====================================== End of "movement of molecules"

   //============================================= Events of molecules
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            //----------------------------Amphiphiles' events
            a_bef = a_arr[1][y][x];  // Amphiphiles out of the membrane decaying
            a_arr[1][y][x] = 0;
            for (k = 0; k < a_bef; k++)
            {
                if (randd() < PAD) pa_arr[0][y][x]++;
                else a_arr[0][y][x]++;
            }
            if (m_arr[y][x] > 0)    // Amphiphiles within the membrane decaying
            {
                m_bef = m_arr[y][x];
                for (k = 0; k < m_bef; k++)
                {
                    if (randd() < PAD * FDW)
                    {
                        m_arr[y][x]--;
                        if (randd() < 0.5) pa_arr[0][y][x]++;  // Amphiphile precursors going into the protocell
                        else							// Amphiphile precursors going out of the protocell
                        {
                            randcase = randl(4);
                            switch (randcase)
                            {
                            case 0:
                                if (m_arr[y][(N + x - 1) % N] == 0) pa_arr[0][y][(N + x - 1) % N]++;
                                else pa_arr[0][y][x]++;
                                break;
                            case 1:
                                if (m_arr[y][(x + 1) % N] == 0) pa_arr[0][y][(x + 1) % N]++;
                                else pa_arr[0][y][x]++;
                                break;
                            case 2:
                                if (m_arr[(N + y - 1) % N][x] == 0) pa_arr[0][(N + y - 1) % N][x]++;
                                else pa_arr[0][y][x]++;
                                break;
                            case 3:
                                if (m_arr[(y + 1) % N][x] == 0) pa_arr[0][(y + 1) % N][x]++;
                                else pa_arr[0][y][x]++;
                                break;
                            default: printf("amphiphile within membrane decaying and pa going out of cell error");
                            }
                        }
                    }
                }
            }

            //----------------------- Peptides' events
            for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
            {
                if (pep->men == 0 && pep->length == 1)
                {
                    if (randd() < PAAD)  // Amino acid decaying
                    {
                        pep->prior->next = pep->next;
                        if (pep->next != room_pep_head[0][y][x])  pep->next->prior = pep->prior;
                        paa_arr[0][y][x]++;
                        pep3 = pep;
                        pep = pep->prior;
                        free(pep3);
                    }
                }
                else if (pep->men == 0 && pep->length == 2)
                {
                    if (randd() < PAADE)  // Peptide's end-decaying
                    {
                        paa_arr[0][y][x]++;
                        pep->information[1] = 0;
                        pep->length--;
                    }
                    if (pep->length > 1)
                    {
                        if (randd() < PPBB)   // Peptide bond breaking
                        {
                            pep1 = (struct peptide*)malloc(LEN_pep);
                            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                            pep1->information[0] = pep->information[1];
                            pep1->length = 1;
                            pep1->men = 0;
                            pep1->loc = -1;
                            pep1->prior = room_pep_head[0][y][x];
                            pep1->next = room_pep_head[0][y][x]->next;
                            if (pep1->next != room_pep_head[0][y][x])(pep1->next)->prior = pep1;
                            room_pep_head[0][y][x]->next = pep1;

                            pep->information[1] = 0;
                            pep->length--;
                        }
                    }
                }
                else if (pep->men == 1 && pep->length == 1)
                {
                    if (randd() < PAAD * FDW)
                    {
                        pep->prior->next = pep->next;
                        if (pep->next != room_pep_head[0][y][x])  pep->next->prior = pep->prior;
                        paa_arr[0][y][x]++;
                        pep3 = pep;
                        pep = pep->prior;
                        free(pep3);
                    }
                }
                else if (pep->men == 1 && pep->length == 2)
                {
                    if (randd() < PAADE * FDW)  // Peptide's end-decaying within membrane
                    {
                        paa_arr[0][y][x]++;
                        pep->men = 0;
                        pep->information[1] = 0;
                        pep->length--;
                    }
                    if (pep->length > 1)
                    {
                        if (randd() < PPBB * FDW)  // Peptide bond breaking within membrane
                        {
                            pep->men = 0;
                            pep1 = (struct peptide*)malloc(LEN_pep);
                            if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                            pep1->information[0] = pep->information[1];
                            pep1->length = 1;
                            pep1->men = 0;
                            pep1->loc = -1;
                            pep1->prior = room_pep_head[0][y][x];
                            pep1->next = room_pep_head[0][y][x]->next;
                            if (pep1->next != room_pep_head[0][y][x])(pep1->next)->prior = pep1;
                            room_pep_head[0][y][x]->next = pep1;

                            pep->information[1] = 0;
                            pep->length--;
                        }
                    }
                }
            }
  
            rna_shuffle(1);
            //------------------------------------- RNA's events
            for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next)
            {
                if (p->rna_pep_head->next != p->rna_pep_head) continue;  // Only considering RNA without binding peptides for the following events
                // Random chain ligating
                for (p3 = p->next; p3 != p; p3 = p3->next)
                {
                    if (p3 == room_head[1][y][x]) { p3 = room_head[1][y][x]->next; if (p3 == p)break; }
                    //new
                    if (p3->rna_pep_head->next != p3->rna_pep_head) continue; // Only considering RNA without binding amino acids or peptides
                    if (p3->length2 == 0)
                    {
                        if (randd() < PRL / (p->length1 + p3->length1)) // Random ligating should be more difficult for longer chains.
                        {
                            if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1) { over_max_len++; continue; }
                            for (a = 0; a < p3->length1; a++) p->information[0][a + p->length1] = p3->information[0][a];
                            p->information[0][p->length1 + p3->length1] = 0;
                            p->length1 = p->length1 + p3->length1;
                            (p3->prior)->next = p3->next;
                            if (p3->next != room_head[1][y][x])(p3->next)->prior = p3->prior;

                            free(p3->rna_pep_head);
                            free(p3);
                            break;
                        }
                    }
                }

                // Decaying and degradating
                if (p->length1 == 1)  // Nucleotide decaying
                {
                    if (m_arr[y][x] == 0) f = PND * FDO;
                    else f = PND;
                    if (p->length2 == 0 && randd() < f)
                    {
                        pn_arr[0][y][x]++;
                        (p->prior)->next = p->next;
                        if (p->next != room_head[1][y][x])(p->next)->prior = p->prior;
                        p3 = p;
                        p = p->prior;

                        free(p3->rna_pep_head);
                        free(p3); continue;
                    }
                }
                else               // RNA Degradating
                {                  // Nucleotide residue decaying at the end of an RNA
                    if (m_arr[y][x] == 0) f = PNDE * FDO;
                    else f = PNDE;
                    if (p->length1 > p->length2 && randd() < f) 
                    {
                        pn_arr[0][y][x]++;
                        p->information[0][p->length1 - 1] = 0;
                        p->length1--;
                    }

                    if (p->length1 != 1)  // Phosphodiester bond breaking
                    {
                        if (m_arr[y][x] == 0) f1 = PBB * FDO;
                        else f1 = PBB;
                        for (j = p->length1; j > 1; j--)
                        {
                            f = f1;
                            if (j <= p->length2 && p->nick[j - 1] == 0) f = f * f;  // Breaking of double chain should be more difficult

                            if (randd() < f)
                            {
                                p3 = (struct rna*)malloc(LEN);
                                if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }

                                for (b = 0; b < p->length1 - j + 1; b++) p3->information[0][b] = p->information[0][b + j - 1];
                                p3->information[0][p->length1 - j + 1] = 0;
                                p->information[0][j - 1] = 0;
                                p3->length1 = p->length1 - j + 1;
                                p->length1 = j - 1;

                                for (i = 0; i < MAX_RNA_LENGTH; i++) p3->nick[i] = 0;

                                pep1 = (struct peptide*)malloc(LEN_pep);
                                if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                                p3->rna_pep_head = pep1;
                                pep1->next = p3->rna_pep_head;

                                if (p->length2 > j - 1)
                                {
                                    for (b = 0; b < p->length2 - j + 1; b++)	p3->information[1][b] = p->information[1][b + j - 1];
                                    p3->information[1][p->length2 - j + 1] = 0;
                                    p->information[1][j - 1] = 0;
                                    p3->length2 = p->length2 - j + 1;
                                    p->length2 = j - 1;
                                }
                                else
                                {
                                    p3->information[1][0] = 0;
                                    p3->length2 = 0;
                                }

                                if (p3->length2 != 0)
                                {
                                    p3->nick[0] = 1;
                                    for (a = 1; a < p3->length2; a++)
                                        p3->nick[a] = p->nick[j + a - 1];
                                }

                                p3->prior = room_head[0][y][x];
                                p3->next = room_head[0][y][x]->next;
                                if (p3->next != room_head[0][y][x])(p3->next)->prior = p3;
                                room_head[0][y][x]->next = p3;
                                break;
                            }
                        }
                    }
                }

                //Template-directed attracting
                if (p->length1 > 3 && (p->length2 != 0 || randd() < 0.5)) // p->length1>3: If an RNA is not longer than 3 nt, it would not serve as a template
                {                                                         // randd()<0.5: A single chain RNA turning to template
                    for (p3 = p->next; p3 != p; p3 = p3->next)   //Template-directed attraction of substrates
                    {
                        if (p3 == room_head[1][y][x]) { p3 = room_head[1][y][x]->next; if (p3 == p)break; }

                        if (p3->rna_pep_head->next != p3->rna_pep_head) continue;
                        if (p3->length2 == 0)
                        {
                            if (p3->length1 <= p->length1 - p->length2)
                            {
                                for (flag1 = 0, b = 0; b < p3->length1; b++)
                                {
                                    if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][p->length2 + b]) == 5)continue;
                                    else if (randd() < PFP)continue;  // False base-pairing
                                    else { flag1 = 1; break; }
                                }
                                if (flag1 == 0)
                                {
                                    rtdaddphili = randd();
                                    if (rtdaddphili < PAT)
                                    {
                                        for (a = 0; a < p3->length1; a++)
                                            p->information[1][p->length2 + a] = p3->information[0][p3->length1 - 1 - a];
                                        p->information[1][p->length2 + p3->length1] = 0;

                                        p->nick[p->length2] = 1;
                                        for (a = 1; a < p3->length1; a++)
                                            p->nick[p->length2 + a] = 0;
                                        p->length2 = p->length2 + p3->length1;

                                        (p3->prior)->next = p3->next;
                                        if (p3->next != room_head[1][y][x])(p3->next)->prior = p3->prior;

                                        free(p3->rna_pep_head);
                                        free(p3);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                //Template-directed ligating
                for (a = 1; a < p->length2; a++)
                {
                    if (p->nick[a] == 0) continue;
                    rtdaddlig = randd();
                    if (rtdaddlig < PTL) p->nick[a] = 0;
                }

                // Double chain separating
                if (p->length2 != 0)
                {
                    a = p->length2 - 1;
                    while (p->nick[a] == 0) a--;

                    if (randd() < pow(PSP, sqrt((p->length2 - a) * 1.0)))     // sqrt: considering single chain folding factor
                    {
                        p3 = (struct rna*)malloc(LEN);
                        if (!p3) { printf("\t%dsep--memeout\n", i); exit(0); }

                        for (b = 0; b < p->length2 - a; b++)
                            p3->information[0][b] = p->information[1][p->length2 - 1 - b];
                        p->information[1][a] = 0;

                        p3->information[0][p->length2 - a] = 0;
                        p3->information[1][0] = 0;

                        p3->length1 = p->length2 - a;
                        p->length2 = a;

                        p3->length2 = 0;

                        for (i = 0; i < MAX_RNA_LENGTH; i++) p3->nick[i] = 0;

                        pep1 = (struct peptide*)malloc(LEN_pep);
                        if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                        p3->rna_pep_head = pep1;
                        pep1->next = p3->rna_pep_head;

                        p3->prior = room_head[1][y][x];
                        p3->next = room_head[1][y][x]->next;
                        if (p3->next != room_head[1][y][x])(p3->next)->prior = p3;
                        room_head[1][y][x]->next = p3;
                    }
                }
            }

            for (p = room_head[1][y][x]->next; p != room_head[1][y][x]; p = p->next) fresh_rna(1);

            //------------------------------ Amphiphile precursors' events
            pa_bef = pa_arr[1][y][x];
            pa_arr[1][y][x] = 0;
            for (k = 0; k < pa_bef; k++)
            {
                ramsyn = randd();
                flagamsyn = 1;
                if (ramsyn < PAF) flagamsyn = 0;
                if (flagamsyn == 0) a_arr[0][y][x]++; // Forming an amphiphile
                else pa_arr[0][y][x]++;
            }

            //--------------------------------  Nucleotide precursors' events
            nt_turn = 0;
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
            {									// Finding nucleotide synthetase ribozymes in the room
                 if (findseq_loc2(nsr_seq, nsrlength, p) != -1 && p->length2 == 0)
                    nt_turn++;        // Caculating possible catalytic turns (rounds) in the room
            }
            pn_bef = pn_arr[1][y][x];
            pn_arr[1][y][x] = 0;
            for (k = 0; k < pn_bef; k++)
            {
                rntsyn = randd();
                flagntsyn = 1;
                if (rntsyn < PNF) flagntsyn = 0; // Not catalyzed
                else if (nt_turn > 0)
                {
                    nt_turn--;
                    if (rntsyn < PNFR) flagntsyn = 0; // Catalyzed by NSR
                }

                if (flagntsyn == 0)   // Nucleotide forming
                {				      
                    p3 = (struct rna*)malloc(LEN);
                    if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(0); }
                    randnt = randl(4) + 1;
                    switch (randnt)
                    {
                    case 1:  p3->information[0][0] = A; break;
                    case 2:  p3->information[0][0] = C; break;
                    case 3:  p3->information[0][0] = G; break;
                    case 4:  p3->information[0][0] = U; break;
                    default: printf("form randnt error");
                    }
                    p3->information[0][1] = 0;
                    p3->information[1][0] = 0;
                    p3->length1 = 1;
                    p3->length2 = 0;
                    p3->prior = room_head[0][y][x];
                    p3->next = room_head[0][y][x]->next;
                    if (p3->next != room_head[0][y][x])(p3->next)->prior = p3;
                    room_head[0][y][x]->next = p3;

                    for (i = 0; i < MAX_RNA_LENGTH; i++) p3->nick[i] = 0;

                    pep1 = (struct peptide*)malloc(LEN_pep);
                    if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                    p3->rna_pep_head = pep1;
                    pep1->next = p3->rna_pep_head;
                }
                else pn_arr[0][y][x]++;
            }

            //--------------------------------  Amino acid precursors' events
            paa_bef = paa_arr[1][y][x];
            paa_arr[1][y][x] = 0;
            for (k = 0; k < paa_bef; k++)
            {
                if (randd() < PAAF)
                {
                    pep1 = (struct peptide*)malloc(LEN_pep);
                    if (!pep1) { printf("\tinit1--memeout\n"); exit(0); }
                    randnt = randl(6) + 1;
                    switch (randnt)
                    {
                    case 1:pep1->information[0] = P; pep1->information[1] = 0; break;
                    case 2:pep1->information[0] = Q; pep1->information[1] = 0; break;
                    case 3:pep1->information[0] = R; pep1->information[1] = 0; break;
                    case 4:pep1->information[0] = S; pep1->information[1] = 0; break;
                    case 5:pep1->information[0] = T; pep1->information[1] = 0; break;
                    case 6:pep1->information[0] = L; pep1->information[1] = 0; break;
                    default: printf("form random amino acids error");
                    }
                    pep1->length = 1;
                    pep1->men = 0;
                    pep1->loc = -1;
                    pep1->prior = room_pep_head[0][y][x];
                    pep1->next = room_pep_head[0][y][x]->next;
                    if (pep1->next != room_pep_head[0][y][x])(pep1->next)->prior = pep1;
                    room_pep_head[0][y][x]->next = pep1;
                }
                else paa_arr[0][y][x]++;
            }
        }
    }
    //======================= End of "events of molecules"

    //======================================== Formation and break of protocells
    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            if (m_arr[y][x] == 0)    // No membrane
            {					 
                if (a_arr[0][y][x] >= LAM)
                {
                    if (randd() < 1 - pow(1 - PMF, a_arr[0][y][x] - LAM + 1))  // Protocell forming
                    {
                        m_arr[y][x] = a_arr[0][y][x];
                        a_arr[0][y][x] = 0;
                    }
                }
            }
            else     // Membrane existing
            {		 
                if (m_arr[y][x] < LAM || randd() < PCB)  // Protocell breaking
                {
                    a_arr[0][y][x] += m_arr[y][x];
                    m_arr[y][x] = 0;
                    for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
                        if (pep->men == 1) pep->men = 0;
                }
            }
        }
    }
    //========================= End of cell forming and breaking

    //==================================== Movement of protocells (including cell fusion)
    for (y = 0; y < N; y++)   // Initially flagging all rooms as "having not been considered yet"
        for (x = 0; x < N; x++)
            flagcmov[y][x] = 0;
    avail_xy_init();      // Initialization for xy_choose
    for (d = 0; d < ROOMNUM; d++)
    {
        xy_choose();    // Picks a room at random
        if (m_arr[y][x] > 0 && flagcmov[y][x] == 0 && randd() < PMC)
        {
            randcase = randl(4);   // Four possible directions
            switch (randcase)
            {
            case 0:  // To left
                if (m_arr[y][(N + x - 1) % N] > 0) // There is a protocell at left
                {						  // Possibly fusing to Left
                    if (randd() < PCF)
                    {
                        cell_fusing(0, -1);
                        flagcmov[y][(N + x - 1) % N] = 1;
                    }
                }
                else    // There is not a protocell at left
                {		// Possibly moving to Left
                    left = 1; up = 1; down = 1;
                    if (m_arr[y][(N + x - 2) % N] > 0) left = 0;
                    if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) up = 0;
                    if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) down = 0;
                    if (left == 0 && up == 0 && down == 0)break; // No way for outer molecules moving, so no cell moving
                    if (left == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, -1, -1);  //  Outer molecules moving allowed only to left
                    if (left == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, -1, 0); // Only up
                    if (left == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, -1, 0); // Only down
                    if (left == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, -1, -1, 0); // Left and up
                    if (left == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, -1, -1, 0); // Left and down
                    if (left == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, -1, 0, 0); // Up and down
                    if (left == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, -1, -1, 0, 0); // Left, up and down
                    cell_moving(0, -1); // Moving to left
                    flagcmov[y][(N + x - 1) % N] = 1; // Flagging this room as "having been considered already"
                }
                break;
            case 1:    // To right
                if (m_arr[y][(x + 1) % N] > 0)
                {
                    if (randd() < PCF)
                    {
                        cell_fusing(0, 1);
                        flagcmov[y][(x + 1) % N] = 1;
                    }
                }
                else
                {
                    right = 1; up = 1; down = 1;
                    if (m_arr[y][(x + 2) % N] > 0) right = 0;
                    if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) up = 0;
                    if (m_arr[(y + 1) % N][(x + 1) % N] > 0) down = 0;
                    if (right == 0 && up == 0 && down == 0)break;
                    if (right == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, 1, 1);
                    if (right == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, 1, 0);
                    if (right == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, 1, 0);
                    if (right == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, 1, 1, 0);
                    if (right == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, 1, 1, 0);
                    if (right == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, 1, 0, 0);
                    if (right == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, 1, 1, 0, 0);
                    cell_moving(0, 1);
                    flagcmov[y][(x + 1) % N] = 1;
                }
                break;
            case 2:    // To up
                if (m_arr[(N + y - 1) % N][x] > 0)
                {
                    if (randd() < PCF)
                    {
                        cell_fusing(-1, 0);
                        flagcmov[(N + y - 1) % N][x] = 1;
                    }
                }
                else
                {
                    up = 1; left = 1; right = 1;
                    if (m_arr[(N + y - 2) % N][x] > 0) up = 0;
                    if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) left = 0;
                    if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) right = 0;
                    if (up == 0 && left == 0 && right == 0)break;
                    if (up == 1 && left == 0 && right == 0) oneway_outer_moving(-1, -1, 0, 0);
                    if (up == 0 && left == 1 && right == 0) oneway_outer_moving(-1, 0, 0, -1);
                    if (up == 0 && left == 0 && right == 1) oneway_outer_moving(-1, 0, 0, 1);
                    if (up == 1 && left == 1 && right == 0) twoway_outer_moving(-1, -1, 0, 0, 0, -1);
                    if (up == 1 && left == 0 && right == 1) twoway_outer_moving(-1, -1, 0, 0, 0, 1);
                    if (up == 0 && left == 1 && right == 1) twoway_outer_moving(-1, 0, 0, 0, -1, 1);
                    if (up == 1 && left == 1 && right == 1) threeway_outer_moving(-1, -1, 0, 0, 0, 0, -1, 1);
                    cell_moving(-1, 0);
                    flagcmov[(N + y - 1) % N][x] = 1;
                }
                break;
            case 3:    // To down
                if (m_arr[(y + 1) % N][x] > 0)
                {
                    if (randd() < PCF)
                    {
                        cell_fusing(1, 0);
                        flagcmov[(y + 1) % N][x] = 1;
                    }
                }
                else
                {
                    down = 1; left = 1; right = 1;
                    if (m_arr[(y + 2) % N][x] > 0) down = 0;
                    if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) left = 0;
                    if (m_arr[(y + 1) % N][(x + 1) % N] > 0) right = 0;
                    if (down == 0 && left == 0 && right == 0)break;
                    if (down == 1 && left == 0 && right == 0) oneway_outer_moving(1, 1, 0, 0);
                    if (down == 0 && left == 1 && right == 0) oneway_outer_moving(1, 0, 0, -1);
                    if (down == 0 && left == 0 && right == 1) oneway_outer_moving(1, 0, 0, 1);
                    if (down == 1 && left == 1 && right == 0) twoway_outer_moving(1, 1, 0, 0, 0, -1);
                    if (down == 1 && left == 0 && right == 1) twoway_outer_moving(1, 1, 0, 0, 0, 1);
                    if (down == 0 && left == 1 && right == 1) twoway_outer_moving(1, 0, 0, 0, -1, 1);
                    if (down == 1 && left == 1 && right == 1) threeway_outer_moving(1, 1, 0, 0, 0, 0, -1, 1);
                    cell_moving(1, 0);
                    flagcmov[(y + 1) % N][x] = 1;
                }
                break;
            default: printf("cell moving error");
            }
        }
    }
    //========================= End of "Movement of protocells (including cell fusion)"

    //========================================== Division of protocells
    for (y = 0; y < N; y++)   // Initially flagging all rooms as "having not been considered yet"
        for (x = 0; x < N; x++)
            flagcdiv[y][x] = 0;
    avail_xy_init();      // Initialization for xy_choose
    for (d = 0; d < ROOMNUM; d++)
    {
        xy_choose();    // Picks a room at random
        if (m_arr[y][x] > 0 && flagcdiv[y][x] == 0 && randd() < PCD * (1 - 2.0 * LAM / m_arr[y][x]))
        {
            randcase = randl(4);   // Four possible directions
            switch (randcase)
            {
            case 0:    // To left
                if (m_arr[y][(N + x - 1) % N] > 0) break;
                left = 1; up = 1; down = 1;
                if (m_arr[y][(N + x - 2) % N] > 0) left = 0;
                if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) up = 0;
                if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) down = 0;
                if (left == 0 && up == 0 && down == 0)break; // No way for outer molecules moving, so no cell division
                if (left == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, -1, -1); // Outer molecules moving allowed only to left
                if (left == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, -1, 0); // Only up
                if (left == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, -1, 0);  // Only down
                if (left == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, -1, -1, 0); // Left and up
                if (left == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, -1, -1, 0);  // Left and down
                if (left == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, -1, 0, 0);  // Up and down
                if (left == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, -1, -1, 0, 0); // Left, up and down
                cell_dividing(0, -1);      // Deviding to left
                flagcdiv[y][(N + x - 1) % N] = 1; // Flagging this room as "having been considered already"
                break;
            case 1: // To right
                if (m_arr[y][(x + 1) % N] > 0) break;
                right = 1; up = 1; down = 1;
                if (m_arr[y][(x + 2) % N] > 0) right = 0;
                if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) up = 0;
                if (m_arr[(y + 1) % N][(x + 1) % N] > 0) down = 0;
                if (right == 0 && up == 0 && down == 0)break;
                if (right == 1 && up == 0 && down == 0) oneway_outer_moving(0, 0, 1, 1);
                if (right == 0 && up == 1 && down == 0) oneway_outer_moving(0, -1, 1, 0);
                if (right == 0 && up == 0 && down == 1) oneway_outer_moving(0, 1, 1, 0);
                if (right == 1 && up == 1 && down == 0) twoway_outer_moving(0, 0, -1, 1, 1, 0);
                if (right == 1 && up == 0 && down == 1) twoway_outer_moving(0, 0, 1, 1, 1, 0);
                if (right == 0 && up == 1 && down == 1) twoway_outer_moving(0, -1, 1, 1, 0, 0);
                if (right == 1 && up == 1 && down == 1) threeway_outer_moving(0, 0, -1, 1, 1, 1, 0, 0);
                cell_dividing(0, 1);
                flagcdiv[y][(x + 1) % N] = 1;
                break;
            case 2: // To up
                if (m_arr[(N + y - 1) % N][x] > 0) break;
                up = 1; left = 1; right = 1;
                if (m_arr[(N + y - 2) % N][x] > 0) up = 0;
                if (m_arr[(N + y - 1) % N][(N + x - 1) % N] > 0) left = 0;
                if (m_arr[(N + y - 1) % N][(x + 1) % N] > 0) right = 0;
                if (up == 0 && left == 0 && right == 0)break;
                if (up == 1 && left == 0 && right == 0) oneway_outer_moving(-1, -1, 0, 0);
                if (up == 0 && left == 1 && right == 0) oneway_outer_moving(-1, 0, 0, -1);
                if (up == 0 && left == 0 && right == 1) oneway_outer_moving(-1, 0, 0, 1);
                if (up == 1 && left == 1 && right == 0) twoway_outer_moving(-1, -1, 0, 0, 0, -1);
                if (up == 1 && left == 0 && right == 1) twoway_outer_moving(-1, -1, 0, 0, 0, 1);
                if (up == 0 && left == 1 && right == 1) twoway_outer_moving(-1, 0, 0, 0, -1, 1);
                if (up == 1 && left == 1 && right == 1) threeway_outer_moving(-1, -1, 0, 0, 0, 0, -1, 1);
                cell_dividing(-1, 0);
                flagcdiv[(N + y - 1) % N][x] = 1;
                break;
            case 3:    // To down
                if (m_arr[(y + 1) % N][x] > 0) break;
                down = 1; left = 1; right = 1;
                if (m_arr[(y + 2) % N][x] > 0) down = 0;
                if (m_arr[(y + 1) % N][(N + x - 1) % N] > 0) left = 0;
                if (m_arr[(y + 1) % N][(x + 1) % N] > 0) right = 0;
                if (down == 0 && left == 0 && right == 0)break;
                if (down == 1 && left == 0 && right == 0) oneway_outer_moving(1, 1, 0, 0);
                if (down == 0 && left == 1 && right == 0) oneway_outer_moving(1, 0, 0, -1);
                if (down == 0 && left == 0 && right == 1) oneway_outer_moving(1, 0, 0, 1);
                if (down == 1 && left == 1 && right == 0) twoway_outer_moving(1, 1, 0, 0, 0, -1);
                if (down == 1 && left == 0 && right == 1) twoway_outer_moving(1, 1, 0, 0, 0, 1);
                if (down == 0 && left == 1 && right == 1) twoway_outer_moving(1, 0, 0, 0, -1, 1);
                if (down == 1 && left == 1 && right == 1) threeway_outer_moving(1, 1, 0, 0, 0, 0, -1, 1);
                cell_dividing(1, 0);
                flagcdiv[(y + 1) % N][x] = 1;
                break;
            default: printf("cell division error");
            }
        }
    }
    //=========================== End of "Division of protocells"
}

/////////////////////////////////////////////////////////////////////////
void record(void)        // Data recording at every interval step (RECINT)
{
    FILE* fptxt;
    errno_t err;
    err = fopen_s(&fptxt, RECTXT, "at");
    if (err != 0) { printf("cannot open file");  exit(-1); }

    total_nt_mat[g] = 0;   // Total materials for nucleotide precursors and RNA (quotients in measurement of nucleotides)
    total_am_mat[g] = 0;   // Total materials for amphiphile precursors and amphiphiles (quotients in measurement of amphiphiles)
    pn_num[g] = 0;         // Number of nucleotide precursors
    rna_num[g] = 0;        // Number of RNA
    pam_num[g] = 0;        // Number of amphiphile precursors
    am_num[g] = 0;         // Number of amphiphiles

    total_pep_num[g] = 0; // Total materials for amino acid precursors and peptides (quotients in measurement of amino acids)
    paa_num[g] = 0;       // Number of amino acid precursors
    pep_num[g] = 0;       // Number of peptides
    STpepM_num[g] = 0;    // Number of MSP within the membrane
    STpepOM_num[g] = 0;   // Number of MSP out of the membrane
    PRpep_num[g] = 0;     // Number of the control peptide

    STpep_rna_num[g] = 0;  // Number of the RNA-gene encoding the MSP
    PRpep_rna_num[g] = 0;  // Number of the RNA-gene encoding the control peptide
    nsr_num[g] = 0;        // Number of NSR
    ctrl_num[g] = 0;       // Number of the control RNA

    cell_num[g] = 0;       // Number of protocells
    cell_nsr_num[g] = 0;   // Number of the protocells containing NSR
    cell_st_num[g] = 0;	   // Number of the protocells containing MSP
    cell_stnsr_num[g] = 0;  // Number of the protocells containing both MSP and NSR
    cell_pr_num[g] = 0;     // Number of the protocells containing the RNA-gene encoding the control peptide
    cell_ctrl_num[g] = 0;   // Number of the protocells containing the control RNA

    for (y = 0; y < N; y++)
    {
        for (x = 0; x < N; x++)
        {
            pn_num[g] += pn_arr[0][y][x];
            pam_num[g] += pa_arr[0][y][x];
            am_num[g] += a_arr[0][y][x];
            total_am_mat[g] += pa_arr[0][y][x] + a_arr[0][y][x] + m_arr[y][x];
            paa_num[g] += paa_arr[0][y][x];
            total_pep_num[g] += paa_arr[0][y][x];
            flagst = 0, flagnsr = 0, flagpr = 0, flagctrl = 0;
            for (pep = room_pep_head[0][y][x]->next; pep != room_pep_head[0][y][x]; pep = pep->next)
            {
                pep_num[g]++;
                total_pep_num[g] += pep->length;
                if (pep->length == 2 && (pep->information[0] == S && pep->information[1] == T))
                {
                    if (pep->men == 1) STpepM_num[g]++;  // MSP within the membrane
                    else STpepOM_num[g]++;     // MSP out of the membrane
                }
                if (pep->length == 2 && (pep->information[0] == P && pep->information[1] == R))
                {
                    PRpep_num[g]++;
                }
            }
            for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
            {
                rna_num[g]++;
                total_nt_mat[g] += p->length1 + p->length2;

                for (pep = p->rna_pep_head->next; pep != p->rna_pep_head; pep = pep->next)
                {
                    pep_num[g]++;
                    total_pep_num[g] += pep->length;
                }

                if (findseq_loc2(st_codseq, st_codlength, p) != -1)
                {
                    STpep_rna_num[g]++;
                    if (flagst == 0) flagst = 1;
                }
                if (findseq_loc2(nsr_seq, nsrlength, p) != -1)
                {
                    nsr_num[g]++;
                    if (flagnsr == 0) flagnsr = 1;
                }
                if (findseq_loc2(pr_codseq, pr_codlength, p) != -1)
                {
                    PRpep_rna_num[g]++;
                    if (flagpr == 0) flagpr = 1;
                }

                if (findseq_loc2(ctrl_seq, ctrllength, p) != -1)
                {
                    ctrl_num[g]++;
                    if (flagctrl == 0) flagctrl = 1;
                }
            }
            if (m_arr[y][x] > 0)
            {
                cell_num[g]++;
                if (flagnsr == 1 && flagst == 0)cell_nsr_num[g]++;
                if (flagnsr == 0 && flagst == 1)cell_st_num[g]++;
                if (flagnsr == 1 && flagst == 1)cell_stnsr_num[g]++;
                if (flagpr == 1)cell_pr_num[g]++;
                if (flagctrl == 1)cell_ctrl_num[g]++;
            }
        }
    }
    total_nt_mat[g] += pn_num[g];

    printf("\nStep=%d: STrna=%d, nsr=%d, PRrna=%d, ctrl=%d, STpep=%d, STpepM=%d, PRpep=%d, C=%d, Cst=%d, Cnsr=%d, Cstnsr=%d, Cpr=%d, Cctrl=%d\n\t [tn=%d,pn=%d,rna=%d,ta=%d,pa=%d,a=%d,tpep=%d,paa=%d,pep=%d]\n",
        step, (int)STpep_rna_num[g], (int)nsr_num[g], (int)PRpep_rna_num[g], (int)ctrl_num[g], (int)STpepM_num[g] + (int)STpepOM_num[g], (int)STpepM_num[g], (int)PRpep_num[g],
        (int)cell_num[g], (int)cell_st_num[g], (int)cell_nsr_num[g], (int)cell_stnsr_num[g], (int)cell_pr_num[g], (int)cell_ctrl_num[g],
        (int)total_nt_mat[g], (int)pn_num[g], (int)rna_num[g], (int)total_am_mat[g], (int)pam_num[g], (int)am_num[g], (int)total_pep_num[g], (int)paa_num[g], (int)pep_num[g]);

    fprintf(fptxt, "Step=%d: STrna=%d, nsr=%d, PRrna=%d, ctrl=%d, STpep=%d, STpepM=%d, PRpep=%d, C=%d, Cst=%d, Cnsr=%d, Cstnsr=%d, Cpr=%d, Cctrl=%d \n",
        step, (int)STpep_rna_num[g], (int)nsr_num[g], (int)PRpep_rna_num[g], (int)ctrl_num[g], (int)STpepM_num[g] + (int)STpepOM_num[g], (int)STpepM_num[g], (int)PRpep_num[g],
        (int)cell_num[g], (int)cell_st_num[g], (int)cell_nsr_num[g], (int)cell_stnsr_num[g], (int)cell_pr_num[g], (int)cell_ctrl_num[g]);
    g++;
    fclose(fptxt);
}

/////////////////////////////////////////////////////////////////////////
void freepool(void)        // Memory releasing
{
    int m;
    for (m = 0; m < 2; m++)
    {
        for (y = 0; y < N; y++)
        {
            for (x = 0; x < N; x++)
            {
                while (1)
                {
                    if (room_head[m][y][x]->next != room_head[m][y][x])
                    {
                        p = room_head[m][y][x]->next;
                        room_head[m][y][x]->next = p->next;
                        while (1)
                        {
                            if (p->rna_pep_head->next != p->rna_pep_head)
                            {
                                pep = p->rna_pep_head->next;
                                p->rna_pep_head->next = pep->next;
                                free(pep);
                            }
                            else break;
                        }
                        free(p->rna_pep_head);
                        free(p);
                    }
                    else break;
                }
                free(room_head[m][y][x]);

                while (1)
                {
                    if (room_pep_head[m][y][x]->next != room_pep_head[m][y][x])
                    {
                        pep = room_pep_head[m][y][x]->next;
                        room_pep_head[m][y][x]->next = pep->next;
                        free(pep);
                    }
                    else break;
                }
                free(room_pep_head[m][y][x]);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////
int main()
{
    inits();            // Initialization of the system

    for (step = 0; step <= STEPNUM; step++)  // The Monte-Carlo loop
    {
        if (step == INOCUSTEP) inoculate();  // Inoculation of RNA-based protocells
        if (step >= STAREC && step % RECINT == 0)    // Data recording at every interval step 
        {
            record();
        }
        unit_action();		// Action of units (molecules and protocells) in the system
    }

    freepool();      // Memory releasing
    return (0);
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of the program











