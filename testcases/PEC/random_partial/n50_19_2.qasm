OPENQASM 2.0;
include "qelib1.inc";
qreg q[56];
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
ccx q[1], q[40], q[42];
s q[11];
h q[39];
ccx q[7], q[36], q[4];
ccx q[21], q[28], q[31];
h q[9];
ccx q[42], q[14], q[28];
s q[44];
t q[18];
ccx q[17], q[40], q[9];
t q[13];
s q[24];
s q[36];
cx q[11], q[24];
t q[29];
s q[48];
ccx q[27], q[12], q[4];
ccx q[6], q[2], q[26];
ccx q[15], q[33], q[21];
s q[20];
cx q[22], q[28];
cx q[41], q[33];
t q[40];
s q[30];
h q[2];
ccx q[38], q[33], q[40];
cx q[24], q[7];
cx q[10], q[27];
ccx q[12], q[33], q[42];
cx q[14], q[27];
h q[29];
ccx q[1], q[10], q[36];
ccx q[7], q[12], q[11];
h q[7];
s q[39];
cx q[38], q[42];
ccx q[28], q[32], q[10];
ccx q[40], q[15], q[2];
ccx q[43], q[16], q[24];
s q[31];
s q[35];
ccx q[46], q[11], q[36];
s q[37];
ccx q[10], q[32], q[38];
h q[35];
ccx q[19], q[11], q[42];
s q[10];
cx q[27], q[15];
s q[45];
h q[19];
s q[29];
h q[12];
cx q[39], q[13];
ccx q[28], q[7], q[20];
ccx q[34], q[30], q[1];
s q[10];
t q[0];
s q[36];
ccx q[32], q[7], q[11];
cx q[4], q[33];
ccx q[37], q[22], q[11];
h q[45];
ccx q[22], q[13], q[43];
t q[29];
cx q[12], q[11];
cx q[45], q[20];
h q[45];
cx q[24], q[22];
cx q[13], q[9];
t q[6];
cx q[16], q[17];
ccx q[10], q[28], q[6];
h q[48];
h q[34];
s q[22];
t q[33];
s q[17];
cx q[25], q[37];
s q[21];
cx q[11], q[12];
ccx q[6], q[30], q[4];
cx q[28], q[2];
cx q[41], q[37];
s q[35];
t q[9];
cx q[26], q[37];
s q[45];
t q[45];
ccx q[10], q[48], q[19];
ccx q[46], q[2], q[44];
s q[6];
s q[37];
s q[7];
ccx q[34], q[42], q[10];
cx q[40], q[18];
ccx q[6], q[38], q[15];
ccx q[36], q[1], q[29];
h q[18];
h q[27];
cx q[8], q[40];
t q[34];
s q[49];
ccx q[44], q[4], q[29];
ccx q[44], q[21], q[40];
ccx q[10], q[16], q[21];
cx q[48], q[3];
ccx q[49], q[6], q[31];
t q[7];
cx q[11], q[1];
ccx q[4], q[11], q[5];
t q[43];
t q[17];
s q[29];
h q[12];
cx q[22], q[23];
t q[13];
h q[10];
cx q[25], q[31];
s q[45];
ccx q[34], q[19], q[3];
ccx q[18], q[27], q[45];
t q[42];
s q[6];
cx q[20], q[18];
t q[22];
h q[0];
t q[19];
s q[8];
s q[35];
cx q[42], q[2];
t q[12];
ccx q[45], q[12], q[31];
ccx q[0], q[18], q[29];
ccx q[30], q[49], q[22];
cx q[39], q[7];
ccx q[6], q[30], q[4];
t q[22];
t q[32];
t q[40];
h q[27];
t q[37];
ccx q[24], q[34], q[15];
t q[31];
h q[3];
cx q[28], q[16];
cx q[6], q[15];
ccx q[6], q[4], q[24];
s q[10];
s q[48];
s q[31];
cx q[0], q[1];
tdg q[1];
tdg q[0];
sdg q[0];
cx q[0], q[1];
z q[2];
y q[4];
sdg q[5];
sdg q[7];
cx q[8], q[7];
tdg q[7];
sdg q[9];
sdg q[11];
cx q[12], q[11];
tdg q[11];
cx q[13], q[14];
sdg q[14];
tdg q[14];
cx q[13], q[14];
tdg q[13];
z q[15];
cx q[17], q[16];
sdg q[16];
cx q[18], q[19];
sdg q[19];
sdg q[18];
cx q[18], q[19];
tdg q[18];
cx q[21], q[20];
sdg q[21];
sdg q[20];
cx q[21], q[20];
tdg q[20];
cx q[22], q[23];
tdg q[23];
cx q[22], q[23];
tdg q[22];
x q[24];
s q[26];
s q[26];
h q[43];
t q[25];
h q[25];
cx q[38], q[33];
cx q[37], q[45];
cx q[33], q[39];
t q[41];
ccx q[47], q[44], q[41];
cx q[38], q[49];
t q[26];
s q[35];
h q[32];
t q[48];
h q[42];
s q[39];
t q[48];
h q[29];
ccx q[44], q[37], q[31];
h q[47];
t q[38];
cx q[30], q[41];
h q[26];
h q[43];
cx q[50], q[16];
cx q[51], q[43];
cx q[52], q[29];
cx q[53], q[32];
cx q[54], q[9];
cx q[55], q[24];
