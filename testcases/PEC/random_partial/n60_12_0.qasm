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
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
ccx q[31], q[4], q[19];
t q[54];
t q[40];
cx q[1], q[13];
ccx q[51], q[16], q[2];
ccx q[3], q[42], q[18];
s q[20];
t q[20];
ccx q[59], q[12], q[9];
ccx q[57], q[13], q[49];
s q[39];
s q[19];
t q[41];
h q[45];
cx q[51], q[21];
cx q[41], q[8];
s q[58];
t q[32];
s q[0];
h q[51];
h q[53];
t q[6];
h q[26];
cx q[37], q[18];
s q[58];
h q[24];
ccx q[7], q[20], q[19];
h q[45];
h q[41];
cx q[48], q[0];
cx q[37], q[51];
s q[6];
h q[0];
s q[15];
s q[57];
s q[17];
ccx q[26], q[27], q[57];
t q[52];
cx q[49], q[32];
cx q[20], q[29];
cx q[55], q[27];
h q[55];
h q[46];
s q[10];
ccx q[0], q[47], q[59];
t q[7];
h q[8];
s q[11];
t q[38];
s q[27];
ccx q[24], q[9], q[42];
ccx q[2], q[57], q[54];
cx q[57], q[6];
s q[31];
s q[55];
ccx q[38], q[0], q[2];
h q[11];
ccx q[32], q[50], q[25];
h q[4];
t q[9];
s q[45];
t q[21];
h q[15];
t q[12];
ccx q[2], q[38], q[36];
t q[26];
cx q[21], q[18];
t q[51];
ccx q[28], q[47], q[10];
cx q[49], q[47];
t q[15];
cx q[1], q[8];
h q[32];
h q[13];
cx q[52], q[2];
ccx q[33], q[39], q[7];
cx q[41], q[27];
ccx q[17], q[8], q[16];
t q[20];
t q[13];
ccx q[5], q[41], q[22];
t q[3];
h q[20];
cx q[24], q[47];
ccx q[36], q[59], q[16];
h q[28];
s q[34];
cx q[55], q[43];
ccx q[17], q[52], q[6];
s q[5];
t q[42];
t q[46];
h q[50];
cx q[47], q[0];
h q[52];
h q[44];
s q[58];
s q[23];
cx q[5], q[13];
h q[39];
s q[19];
cx q[17], q[11];
t q[0];
s q[41];
s q[2];
h q[49];
h q[4];
t q[56];
h q[32];
t q[7];
t q[38];
h q[46];
cx q[45], q[40];
t q[5];
s q[44];
ccx q[46], q[19], q[35];
h q[27];
cx q[20], q[55];
cx q[28], q[14];
ccx q[28], q[40], q[27];
ccx q[57], q[13], q[28];
s q[36];
h q[28];
h q[23];
h q[54];
t q[50];
t q[24];
cx q[5], q[44];
cx q[15], q[21];
h q[22];
cx q[12], q[26];
t q[50];
ccx q[38], q[42], q[13];
cx q[49], q[30];
t q[42];
h q[17];
s q[7];
s q[59];
h q[22];
t q[57];
cx q[42], q[3];
cx q[51], q[46];
h q[13];
t q[39];
t q[40];
t q[21];
ccx q[28], q[36], q[44];
s q[34];
s q[1];
t q[39];
ccx q[47], q[40], q[56];
h q[31];
cx q[29], q[45];
cx q[45], q[49];
cx q[25], q[33];
h q[58];
s q[23];
ccx q[44], q[4], q[56];
s q[30];
cx q[57], q[55];
t q[14];
s q[36];
s q[41];
h q[15];
t q[4];
s q[40];
cx q[25], q[28];
s q[54];
s q[50];
ccx q[35], q[10], q[43];
cx q[0], q[52];
s q[43];
ccx q[22], q[46], q[30];
ccx q[27], q[19], q[34];
ccx q[39], q[53], q[49];
cx q[24], q[12];
t q[50];
h q[56];
cx q[9], q[17];
cx q[57], q[36];
cx q[0], q[1];
tdg q[0];
tdg q[1];
cx q[0], q[1];
tdg q[0];
tdg q[3];
tdg q[2];
cx q[2], q[3];
sdg q[2];
sdg q[4];
cx q[4], q[5];
tdg q[5];
cx q[4], q[5];
sdg q[4];
cx q[6], q[7];
tdg q[7];
cx q[6], q[7];
sdg q[6];
cx q[9], q[8];
tdg q[8];
cx q[9], q[8];
tdg q[8];
sdg q[11];
cx q[11], q[10];
tdg q[10];
sdg q[12];
cx q[12], q[13];
cx q[14], q[15];
sdg q[14];
tdg q[17];
cx q[16], q[17];
sdg q[16];
cx q[18], q[19];
sdg q[19];
cx q[18], q[19];
tdg q[18];
sdg q[18];
cx q[21], q[20];
tdg q[20];
cx q[21], q[20];
tdg q[20];
tdg q[23];
sdg q[22];
cx q[22], q[23];
tdg q[24];
cx q[25], q[24];
tdg q[24];
sdg q[26];
cx q[29], q[28];
tdg q[28];
cx q[37], q[52];
t q[33];
cx q[50], q[31];
s q[45];
h q[40];
s q[54];
ccx q[40], q[37], q[33];
t q[46];
cx q[54], q[44];
s q[54];
s q[37];
s q[40];
t q[34];
h q[33];
s q[38];
h q[57];
ccx q[37], q[49], q[48];
s q[52];
t q[54];
h q[31];
cx q[37], q[44];
h q[57];
cx q[37], q[35];
s q[43];
t q[34];
ccx q[48], q[37], q[44];
h q[36];
cx q[42], q[56];
ccx q[49], q[37], q[56];
h q[30];
