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
s q[15];
t q[31];
ccx q[47], q[48], q[30];
ccx q[29], q[27], q[14];
h q[19];
cx q[6], q[3];
ccx q[2], q[33], q[18];
cx q[11], q[19];
ccx q[22], q[38], q[36];
h q[5];
ccx q[10], q[2], q[47];
t q[40];
s q[11];
s q[42];
s q[49];
cx q[8], q[0];
t q[16];
h q[37];
h q[7];
ccx q[1], q[8], q[19];
ccx q[6], q[43], q[19];
cx q[46], q[7];
ccx q[34], q[31], q[39];
h q[43];
h q[34];
h q[35];
s q[15];
ccx q[37], q[39], q[28];
ccx q[35], q[42], q[16];
cx q[4], q[30];
ccx q[18], q[15], q[36];
t q[45];
ccx q[25], q[6], q[27];
h q[26];
ccx q[37], q[4], q[17];
h q[35];
ccx q[6], q[32], q[11];
ccx q[25], q[12], q[33];
t q[13];
t q[22];
t q[16];
t q[29];
h q[36];
t q[13];
ccx q[2], q[13], q[8];
h q[0];
s q[22];
ccx q[11], q[34], q[3];
t q[19];
h q[29];
s q[32];
h q[3];
cx q[3], q[22];
ccx q[7], q[13], q[4];
s q[12];
h q[38];
t q[45];
cx q[46], q[22];
ccx q[49], q[19], q[23];
cx q[12], q[23];
cx q[11], q[10];
s q[49];
cx q[16], q[30];
ccx q[4], q[36], q[19];
t q[33];
h q[7];
ccx q[29], q[31], q[8];
s q[32];
s q[28];
h q[28];
h q[48];
ccx q[26], q[22], q[8];
s q[23];
h q[37];
s q[30];
cx q[30], q[46];
h q[43];
ccx q[26], q[16], q[21];
h q[39];
cx q[5], q[15];
t q[44];
t q[14];
t q[2];
cx q[16], q[19];
ccx q[48], q[28], q[10];
s q[8];
t q[49];
ccx q[3], q[46], q[28];
ccx q[40], q[25], q[46];
ccx q[5], q[30], q[34];
ccx q[2], q[23], q[22];
h q[34];
t q[0];
t q[40];
ccx q[47], q[15], q[10];
s q[4];
cx q[27], q[38];
s q[48];
cx q[43], q[33];
t q[22];
h q[15];
s q[30];
cx q[0], q[41];
cx q[45], q[10];
t q[40];
ccx q[45], q[33], q[38];
cx q[34], q[11];
ccx q[30], q[38], q[37];
ccx q[0], q[6], q[38];
cx q[17], q[12];
cx q[11], q[28];
t q[40];
cx q[44], q[41];
ccx q[35], q[0], q[32];
s q[12];
cx q[16], q[27];
h q[23];
cx q[26], q[32];
ccx q[38], q[1], q[33];
cx q[27], q[31];
ccx q[38], q[39], q[40];
ccx q[6], q[13], q[38];
ccx q[15], q[22], q[43];
cx q[47], q[1];
cx q[44], q[13];
ccx q[37], q[9], q[19];
cx q[40], q[14];
cx q[17], q[26];
s q[38];
h q[41];
h q[5];
ccx q[34], q[21], q[18];
h q[35];
s q[33];
cx q[21], q[9];
s q[19];
t q[14];
ccx q[4], q[6], q[8];
ccx q[20], q[34], q[47];
t q[38];
cx q[21], q[30];
s q[34];
cx q[42], q[7];
t q[39];
ccx q[47], q[43], q[9];
h q[12];
t q[5];
t q[7];
cx q[45], q[36];
cx q[17], q[20];
y q[0];
sdg q[1];
z q[4];
sdg q[5];
cx q[6], q[5];
sdg q[5];
cx q[7], q[8];
cx q[9], q[10];
tdg q[9];
cx q[12], q[11];
sdg q[11];
cx q[12], q[11];
x q[11];
cx q[13], q[14];
sdg q[13];
sdg q[15];
cx q[16], q[15];
sdg q[15];
cx q[16], q[15];
tdg q[15];
cx q[17], q[18];
sdg q[18];
cx q[17], q[18];
x q[17];
x q[19];
tdg q[21];
cx q[20], q[21];
sdg q[20];
cx q[22], q[23];
sdg q[22];
h q[41];
h q[40];
ccx q[43], q[25], q[44];
ccx q[32], q[44], q[26];
cx q[32], q[48];
h q[27];
cx q[33], q[48];
h q[25];
s q[33];
s q[38];
s q[47];
s q[31];
h q[31];
h q[38];
cx q[37], q[27];
h q[42];
ccx q[29], q[31], q[40];
cx q[37], q[32];
ccx q[37], q[33], q[36];
ccx q[44], q[38], q[25];
s q[40];
cx q[39], q[42];
t q[36];
ccx q[27], q[41], q[25];
cx q[28], q[42];
