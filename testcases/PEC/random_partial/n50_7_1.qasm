OPENQASM 2.0;
include "qelib1.inc";
qreg q[50];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[14];
cx q[24], q[46];
s q[30];
h q[40];
t q[5];
t q[42];
h q[17];
cx q[17], q[36];
tdg q[36];
cx q[25], q[36];
t q[36];
cx q[17], q[36];
tdg q[36];
cx q[25], q[36];
t q[36];
cx q[25], q[17];
tdg q[17];
cx q[25], q[17];
t q[25];
t q[17];
h q[17];
h q[32];
cx q[32], q[17];
tdg q[17];
cx q[47], q[17];
t q[17];
cx q[32], q[17];
tdg q[17];
cx q[47], q[17];
t q[17];
cx q[47], q[32];
tdg q[32];
cx q[47], q[32];
t q[47];
t q[32];
h q[32];
cx q[39], q[27];
t q[12];
s q[42];
t q[41];
cx q[7], q[14];
cx q[33], q[9];
h q[13];
cx q[13], q[42];
tdg q[42];
cx q[8], q[42];
t q[42];
cx q[13], q[42];
tdg q[42];
cx q[8], q[42];
t q[42];
cx q[8], q[13];
tdg q[13];
cx q[8], q[13];
t q[8];
t q[13];
h q[13];
h q[31];
cx q[31], q[34];
tdg q[34];
cx q[10], q[34];
t q[34];
cx q[31], q[34];
tdg q[34];
cx q[10], q[34];
t q[34];
cx q[10], q[31];
tdg q[31];
cx q[10], q[31];
t q[10];
t q[31];
h q[31];
h q[20];
h q[41];
cx q[45], q[29];
t q[48];
t q[2];
cx q[39], q[46];
h q[14];
t q[7];
t q[38];
s q[7];
h q[38];
cx q[38], q[6];
tdg q[6];
cx q[7], q[6];
t q[6];
cx q[38], q[6];
tdg q[6];
cx q[7], q[6];
t q[6];
cx q[7], q[38];
tdg q[38];
cx q[7], q[38];
t q[7];
t q[38];
h q[38];
h q[43];
t q[14];
t q[31];
h q[9];
t q[15];
h q[19];
h q[43];
h q[14];
s q[40];
h q[17];
cx q[42], q[26];
h q[48];
cx q[48], q[36];
tdg q[36];
cx q[20], q[36];
t q[36];
cx q[48], q[36];
tdg q[36];
cx q[20], q[36];
t q[36];
cx q[20], q[48];
tdg q[48];
cx q[20], q[48];
t q[20];
t q[48];
h q[48];
cx q[23], q[4];
s q[19];
s q[2];
cx q[6], q[3];
h q[9];
s q[18];
s q[14];
h q[30];
s q[48];
s q[31];
t q[17];
cx q[46], q[37];
cx q[15], q[35];
h q[2];
cx q[2], q[36];
tdg q[36];
cx q[4], q[36];
t q[36];
cx q[2], q[36];
tdg q[36];
cx q[4], q[36];
t q[36];
cx q[4], q[2];
tdg q[2];
cx q[4], q[2];
t q[4];
t q[2];
h q[2];
h q[44];
h q[10];
cx q[10], q[3];
tdg q[3];
cx q[38], q[3];
t q[3];
cx q[10], q[3];
tdg q[3];
cx q[38], q[3];
t q[3];
cx q[38], q[10];
tdg q[10];
cx q[38], q[10];
t q[38];
t q[10];
h q[10];
t q[11];
t q[18];
cx q[25], q[17];
t q[6];
h q[45];
cx q[45], q[46];
tdg q[46];
cx q[30], q[46];
t q[46];
cx q[45], q[46];
tdg q[46];
cx q[30], q[46];
t q[46];
cx q[30], q[45];
tdg q[45];
cx q[30], q[45];
t q[30];
t q[45];
h q[45];
h q[27];
cx q[27], q[33];
tdg q[33];
cx q[42], q[33];
t q[33];
cx q[27], q[33];
tdg q[33];
cx q[42], q[33];
t q[33];
cx q[42], q[27];
tdg q[27];
cx q[42], q[27];
t q[42];
t q[27];
h q[27];
s q[48];
h q[0];
s q[45];
h q[2];
cx q[2], q[20];
tdg q[20];
cx q[45], q[20];
t q[20];
cx q[2], q[20];
tdg q[20];
cx q[45], q[20];
t q[20];
cx q[45], q[2];
tdg q[2];
cx q[45], q[2];
t q[45];
t q[2];
h q[2];
s q[11];
cx q[45], q[41];
h q[48];
cx q[48], q[36];
tdg q[36];
cx q[5], q[36];
t q[36];
cx q[48], q[36];
tdg q[36];
cx q[5], q[36];
t q[36];
cx q[5], q[48];
tdg q[48];
cx q[5], q[48];
t q[5];
t q[48];
h q[48];
t q[39];
cx q[12], q[18];
h q[31];
cx q[31], q[36];
tdg q[36];
cx q[26], q[36];
t q[36];
cx q[31], q[36];
tdg q[36];
cx q[26], q[36];
t q[36];
cx q[26], q[31];
tdg q[31];
cx q[26], q[31];
t q[26];
t q[31];
h q[31];
s q[20];
cx q[26], q[47];
t q[4];
h q[45];
cx q[45], q[29];
tdg q[29];
cx q[41], q[29];
t q[29];
cx q[45], q[29];
tdg q[29];
cx q[41], q[29];
t q[29];
cx q[41], q[45];
tdg q[45];
cx q[41], q[45];
t q[41];
t q[45];
h q[45];
h q[42];
t q[23];
h q[2];
t q[16];
h q[18];
h q[29];
s q[32];
cx q[42], q[41];
t q[16];
cx q[15], q[47];
t q[1];
cx q[7], q[22];
h q[24];
t q[37];
cx q[21], q[49];
cx q[9], q[31];
s q[0];
h q[25];
cx q[11], q[15];
cx q[8], q[41];
h q[7];
h q[46];
t q[38];
t q[28];
cx q[22], q[8];
cx q[11], q[4];
cx q[0], q[12];
h q[9];
s q[14];
h q[23];
s q[14];
t q[45];
t q[37];
s q[33];
t q[28];
s q[21];
h q[44];
s q[21];
h q[20];
cx q[20], q[41];
tdg q[41];
cx q[19], q[41];
t q[41];
cx q[20], q[41];
tdg q[41];
cx q[19], q[41];
t q[41];
cx q[19], q[20];
tdg q[20];
cx q[19], q[20];
t q[19];
t q[20];
h q[20];
s q[8];
h q[2];
cx q[2], q[3];
tdg q[3];
cx q[13], q[3];
t q[3];
cx q[2], q[3];
tdg q[3];
cx q[13], q[3];
t q[3];
cx q[13], q[2];
tdg q[2];
cx q[13], q[2];
t q[13];
t q[2];
h q[2];
s q[42];
cx q[41], q[14];
s q[18];
s q[3];
h q[15];
cx q[15], q[41];
tdg q[41];
cx q[30], q[41];
t q[41];
cx q[15], q[41];
tdg q[41];
cx q[30], q[41];
t q[41];
cx q[30], q[15];
tdg q[15];
cx q[30], q[15];
t q[30];
t q[15];
h q[15];
s q[43];
h q[37];
cx q[37], q[35];
tdg q[35];
cx q[6], q[35];
t q[35];
cx q[37], q[35];
tdg q[35];
cx q[6], q[35];
t q[35];
cx q[6], q[37];
tdg q[37];
cx q[6], q[37];
t q[6];
t q[37];
h q[37];
s q[48];
t q[46];
h q[20];
cx q[20], q[31];
tdg q[31];
cx q[27], q[31];
t q[31];
cx q[20], q[31];
tdg q[31];
cx q[27], q[31];
t q[31];
cx q[27], q[20];
tdg q[20];
cx q[27], q[20];
t q[27];
t q[20];
h q[20];
h q[37];
s q[7];
h q[12];
h q[26];
cx q[26], q[30];
tdg q[30];
cx q[11], q[30];
t q[30];
cx q[26], q[30];
tdg q[30];
cx q[11], q[30];
t q[30];
cx q[11], q[26];
tdg q[26];
cx q[11], q[26];
t q[11];
t q[26];
h q[26];
t q[3];
h q[36];
cx q[36], q[45];
tdg q[45];
cx q[4], q[45];
t q[45];
cx q[36], q[45];
tdg q[45];
cx q[4], q[45];
t q[45];
cx q[4], q[36];
tdg q[36];
cx q[4], q[36];
t q[4];
t q[36];
h q[36];
h q[25];
s q[18];
h q[48];
s q[22];
h q[28];
t q[23];
t q[46];
h q[1];
cx q[1], q[20];
tdg q[20];
cx q[15], q[20];
t q[20];
cx q[1], q[20];
tdg q[20];
cx q[15], q[20];
t q[20];
cx q[15], q[1];
tdg q[1];
cx q[15], q[1];
t q[15];
t q[1];
h q[1];
cx q[9], q[11];
h q[14];
cx q[14], q[31];
tdg q[31];
cx q[17], q[31];
t q[31];
cx q[14], q[31];
tdg q[31];
cx q[17], q[31];
t q[31];
cx q[17], q[14];
tdg q[14];
cx q[17], q[14];
t q[17];
t q[14];
h q[14];
h q[22];
h q[43];
h q[5];
h q[15];
t q[21];
t q[49];
cx q[48], q[22];
h q[18];
cx q[18], q[34];
tdg q[34];
cx q[45], q[34];
t q[34];
cx q[18], q[34];
tdg q[34];
cx q[45], q[34];
t q[34];
cx q[45], q[18];
tdg q[18];
cx q[45], q[18];
t q[45];
t q[18];
h q[18];
y q[0];
y q[1];
t q[3];
cx q[2], q[3];
s q[4];
t q[5];
cx q[4], q[5];
s q[7];
cx q[6], q[7];
t q[6];
s q[8];
t q[9];
cx q[9], q[8];
s q[12];
t q[12];
cx q[13], q[12];
z q[14];
s q[17];
cx q[17], q[18];
s q[18];
cx q[17], q[18];
s q[21];
t q[22];
cx q[21], q[22];
x q[23];
t q[30];
s q[34];
t q[46];
h q[35];
s q[49];
t q[43];
h q[49];
cx q[49], q[42];
tdg q[42];
cx q[29], q[42];
t q[42];
cx q[49], q[42];
tdg q[42];
cx q[29], q[42];
t q[42];
cx q[29], q[49];
tdg q[49];
cx q[29], q[49];
t q[29];
t q[49];
h q[49];
t q[26];
h q[48];
cx q[48], q[33];
tdg q[33];
cx q[49], q[33];
t q[33];
cx q[48], q[33];
tdg q[33];
cx q[49], q[33];
t q[33];
cx q[49], q[48];
tdg q[48];
cx q[49], q[48];
t q[49];
t q[48];
h q[48];
cx q[49], q[26];
h q[33];
t q[43];
s q[34];
h q[26];
h q[26];
cx q[26], q[46];
tdg q[46];
cx q[31], q[46];
t q[46];
cx q[26], q[46];
tdg q[46];
cx q[31], q[46];
t q[46];
cx q[31], q[26];
tdg q[26];
cx q[31], q[26];
t q[31];
t q[26];
h q[26];
cx q[29], q[27];
t q[36];
h q[30];
cx q[30], q[26];
tdg q[26];
cx q[32], q[26];
t q[26];
cx q[30], q[26];
tdg q[26];
cx q[32], q[26];
t q[26];
cx q[32], q[30];
tdg q[30];
cx q[32], q[30];
t q[32];
t q[30];
h q[30];
t q[37];
s q[26];
cx q[43], q[37];
t q[26];
h q[27];
h q[49];
cx q[49], q[44];
tdg q[44];
cx q[37], q[44];
t q[44];
cx q[49], q[44];
tdg q[44];
cx q[37], q[44];
t q[44];
cx q[37], q[49];
tdg q[49];
cx q[37], q[49];
t q[37];
t q[49];
h q[49];
t q[28];