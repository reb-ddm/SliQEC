OPENQASM 2.0;
include "qelib1.inc";
qreg q[59];
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
h q[35];
ccx q[42], q[30], q[15];
cx q[24], q[18];
t q[29];
s q[31];
t q[41];
h q[23];
cx q[9], q[29];
s q[21];
t q[17];
cx q[47], q[7];
cx q[8], q[14];
s q[28];
t q[22];
t q[43];
cx q[49], q[26];
cx q[15], q[28];
cx q[20], q[0];
ccx q[19], q[28], q[7];
t q[24];
t q[25];
ccx q[26], q[37], q[9];
h q[1];
t q[1];
cx q[36], q[5];
s q[42];
t q[39];
s q[49];
h q[25];
cx q[16], q[20];
ccx q[1], q[26], q[47];
h q[7];
cx q[47], q[12];
t q[38];
s q[42];
t q[22];
s q[18];
cx q[31], q[40];
cx q[19], q[1];
s q[3];
t q[2];
t q[39];
s q[39];
cx q[14], q[1];
t q[5];
h q[20];
h q[12];
t q[48];
t q[18];
t q[40];
ccx q[35], q[21], q[48];
h q[15];
s q[22];
h q[30];
t q[33];
ccx q[28], q[26], q[9];
t q[8];
s q[34];
cx q[1], q[2];
h q[35];
h q[27];
cx q[16], q[24];
ccx q[18], q[42], q[40];
cx q[7], q[29];
h q[1];
t q[34];
s q[24];
s q[22];
cx q[47], q[38];
s q[46];
ccx q[4], q[11], q[35];
h q[4];
t q[28];
h q[19];
t q[13];
t q[24];
h q[42];
s q[23];
ccx q[48], q[37], q[34];
t q[13];
ccx q[1], q[8], q[7];
cx q[4], q[8];
cx q[11], q[13];
s q[12];
t q[43];
cx q[18], q[7];
h q[8];
t q[47];
s q[41];
ccx q[29], q[33], q[25];
s q[1];
t q[22];
h q[16];
t q[36];
ccx q[27], q[0], q[24];
s q[24];
ccx q[15], q[44], q[10];
h q[9];
s q[4];
t q[47];
ccx q[46], q[4], q[28];
h q[8];
ccx q[38], q[44], q[8];
s q[14];
h q[46];
ccx q[27], q[10], q[41];
h q[47];
s q[9];
ccx q[17], q[22], q[18];
ccx q[32], q[41], q[46];
t q[21];
s q[6];
cx q[8], q[18];
h q[24];
h q[2];
h q[49];
cx q[29], q[19];
h q[28];
h q[45];
h q[32];
ccx q[5], q[40], q[22];
h q[12];
ccx q[9], q[28], q[48];
h q[47];
ccx q[10], q[8], q[32];
cx q[33], q[14];
h q[6];
h q[42];
h q[11];
cx q[28], q[27];
cx q[44], q[6];
h q[38];
h q[39];
cx q[7], q[18];
ccx q[41], q[19], q[1];
t q[13];
h q[9];
cx q[47], q[12];
cx q[44], q[29];
cx q[11], q[12];
s q[48];
t q[29];
h q[26];
t q[20];
t q[29];
h q[21];
s q[6];
cx q[28], q[9];
s q[3];
cx q[46], q[20];
sdg q[0];
sdg q[1];
cx q[1], q[0];
tdg q[0];
tdg q[2];
cx q[3], q[2];
tdg q[2];
sdg q[2];
sdg q[4];
y q[6];
x q[7];
cx q[9], q[8];
tdg q[8];
cx q[10], q[11];
x q[12];
cx q[13], q[14];
sdg q[14];
tdg q[14];
sdg q[13];
cx q[13], q[14];
y q[15];
tdg q[16];
cx q[19], q[18];
tdg q[19];
tdg q[18];
sdg q[21];
tdg q[20];
cx q[21], q[20];
sdg q[20];
tdg q[22];
x q[22];
cx q[23], q[22];
x q[22];
t q[28];
ccx q[38], q[32], q[33];
t q[39];
cx q[33], q[26];
t q[35];
t q[41];
t q[37];
ccx q[34], q[38], q[39];
cx q[28], q[32];
ccx q[44], q[45], q[40];
h q[27];
ccx q[37], q[47], q[41];
s q[48];
t q[39];
t q[45];
t q[26];
ccx q[42], q[29], q[31];
h q[43];
t q[25];
s q[42];
cx q[32], q[29];
s q[36];
t q[46];
cx q[33], q[40];
t q[49];
cx q[50], q[45];
cx q[51], q[36];
cx q[52], q[35];
cx q[53], q[30];
cx q[54], q[12];
cx q[55], q[48];
cx q[56], q[47];
cx q[57], q[24];
cx q[58], q[16];
