OPENQASM 2.0;
include "qelib1.inc";
qreg q[60];
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
h q[27];
s q[43];
cx q[45], q[7];
s q[33];
h q[22];
s q[37];
h q[32];
s q[20];
cx q[5], q[1];
s q[19];
s q[20];
t q[23];
s q[2];
cx q[11], q[15];
s q[28];
cx q[32], q[49];
s q[37];
ccx q[2], q[44], q[10];
h q[34];
h q[24];
cx q[23], q[5];
h q[6];
s q[39];
ccx q[3], q[10], q[14];
t q[31];
ccx q[47], q[17], q[26];
t q[20];
ccx q[22], q[28], q[27];
t q[34];
s q[29];
s q[0];
t q[15];
s q[23];
ccx q[34], q[44], q[5];
h q[17];
s q[44];
h q[24];
t q[17];
h q[23];
h q[35];
t q[25];
cx q[44], q[34];
cx q[4], q[46];
ccx q[36], q[9], q[18];
ccx q[21], q[5], q[30];
s q[19];
h q[6];
ccx q[10], q[4], q[45];
h q[19];
h q[6];
cx q[10], q[32];
t q[29];
t q[45];
h q[48];
cx q[0], q[4];
ccx q[6], q[2], q[40];
s q[25];
s q[40];
t q[13];
s q[10];
cx q[21], q[7];
ccx q[13], q[48], q[45];
cx q[25], q[21];
h q[23];
ccx q[41], q[42], q[20];
t q[7];
s q[44];
t q[7];
h q[27];
h q[31];
cx q[30], q[40];
ccx q[6], q[38], q[2];
cx q[10], q[31];
t q[15];
t q[26];
s q[13];
cx q[20], q[22];
ccx q[27], q[5], q[8];
h q[48];
t q[4];
cx q[38], q[35];
s q[37];
cx q[34], q[7];
ccx q[11], q[33], q[5];
cx q[37], q[41];
s q[30];
t q[7];
h q[23];
s q[33];
h q[37];
t q[44];
s q[36];
cx q[38], q[37];
cx q[27], q[45];
ccx q[25], q[13], q[11];
s q[1];
ccx q[1], q[6], q[2];
s q[27];
s q[47];
t q[29];
s q[36];
s q[25];
t q[37];
cx q[42], q[32];
cx q[15], q[30];
s q[16];
h q[4];
cx q[30], q[11];
s q[12];
cx q[17], q[36];
cx q[10], q[43];
t q[8];
h q[13];
s q[37];
ccx q[14], q[13], q[1];
s q[3];
cx q[20], q[21];
ccx q[18], q[29], q[34];
h q[34];
h q[23];
s q[27];
ccx q[16], q[10], q[40];
t q[23];
ccx q[45], q[9], q[8];
s q[45];
h q[13];
s q[48];
t q[4];
h q[38];
t q[34];
s q[49];
h q[46];
ccx q[5], q[17], q[27];
h q[23];
s q[48];
s q[21];
s q[19];
s q[48];
t q[40];
s q[25];
t q[8];
s q[13];
s q[10];
s q[19];
t q[39];
h q[37];
s q[13];
s q[20];
s q[35];
h q[17];
cx q[0], q[1];
tdg q[0];
cx q[3], q[2];
tdg q[2];
x q[5];
tdg q[7];
sdg q[9];
tdg q[10];
cx q[10], q[9];
sdg q[9];
sdg q[11];
cx q[12], q[11];
tdg q[11];
x q[13];
tdg q[14];
tdg q[16];
y q[18];
x q[19];
cx q[21], q[20];
tdg q[20];
tdg q[23];
cx q[22], q[23];
tdg q[22];
y q[24];
ccx q[49], q[38], q[44];
s q[47];
s q[27];
t q[48];
s q[31];
h q[40];
s q[30];
t q[25];
cx q[36], q[38];
ccx q[39], q[26], q[45];
t q[27];
h q[45];
s q[45];
t q[40];
cx q[34], q[37];
h q[31];
ccx q[28], q[45], q[39];
s q[26];
ccx q[43], q[36], q[39];
ccx q[48], q[39], q[43];
t q[26];
h q[28];
h q[28];
s q[41];
ccx q[47], q[27], q[46];
cx q[50], q[23];
cx q[51], q[26];
cx q[52], q[25];
cx q[53], q[40];
cx q[54], q[38];
cx q[55], q[12];
cx q[56], q[48];
cx q[57], q[2];
cx q[58], q[9];
cx q[59], q[44];
