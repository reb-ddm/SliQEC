OPENQASM 2.0;
include "qelib1.inc";
qreg q[78];
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
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
ccx q[35], q[16], q[65];
ccx q[57], q[44], q[30];
cx q[15], q[49];
s q[44];
s q[18];
h q[1];
h q[39];
cx q[64], q[20];
ccx q[24], q[17], q[28];
h q[63];
t q[36];
cx q[19], q[9];
s q[67];
s q[30];
t q[57];
cx q[69], q[35];
t q[1];
t q[4];
ccx q[69], q[14], q[2];
s q[66];
h q[18];
cx q[35], q[28];
t q[21];
h q[45];
h q[22];
cx q[20], q[38];
ccx q[37], q[0], q[61];
s q[18];
ccx q[7], q[16], q[65];
t q[34];
h q[7];
cx q[46], q[49];
cx q[38], q[59];
h q[6];
t q[38];
s q[10];
ccx q[9], q[53], q[3];
ccx q[65], q[15], q[57];
cx q[41], q[44];
cx q[38], q[64];
cx q[10], q[53];
s q[33];
h q[66];
s q[33];
cx q[6], q[5];
h q[0];
t q[0];
s q[28];
h q[28];
ccx q[48], q[1], q[42];
ccx q[65], q[26], q[9];
h q[44];
s q[13];
ccx q[67], q[68], q[46];
cx q[0], q[18];
ccx q[40], q[50], q[51];
cx q[33], q[40];
h q[26];
h q[16];
cx q[42], q[54];
t q[36];
h q[32];
s q[7];
t q[5];
h q[32];
h q[49];
cx q[19], q[49];
s q[44];
s q[4];
h q[67];
t q[3];
t q[17];
h q[3];
s q[15];
s q[8];
t q[39];
ccx q[23], q[22], q[30];
cx q[13], q[19];
h q[37];
h q[37];
ccx q[6], q[7], q[15];
s q[49];
ccx q[2], q[28], q[24];
t q[44];
cx q[52], q[6];
t q[59];
cx q[25], q[64];
s q[29];
s q[22];
t q[49];
t q[69];
cx q[5], q[61];
ccx q[18], q[4], q[6];
cx q[47], q[57];
s q[4];
ccx q[36], q[20], q[3];
s q[63];
cx q[22], q[66];
t q[54];
h q[25];
s q[48];
s q[57];
s q[47];
h q[65];
h q[10];
ccx q[24], q[7], q[35];
h q[20];
h q[12];
s q[47];
h q[56];
s q[32];
h q[7];
h q[16];
h q[1];
ccx q[1], q[14], q[46];
cx q[59], q[63];
cx q[34], q[22];
s q[23];
s q[37];
h q[40];
t q[51];
t q[35];
s q[2];
t q[17];
s q[48];
h q[47];
cx q[28], q[48];
ccx q[36], q[11], q[8];
t q[17];
s q[43];
ccx q[34], q[68], q[21];
t q[38];
t q[50];
s q[3];
t q[18];
s q[22];
h q[54];
h q[5];
s q[62];
s q[1];
ccx q[9], q[60], q[61];
ccx q[51], q[66], q[13];
h q[59];
h q[12];
ccx q[45], q[56], q[10];
ccx q[10], q[55], q[6];
h q[62];
cx q[31], q[4];
s q[32];
cx q[33], q[69];
ccx q[27], q[13], q[1];
h q[54];
h q[7];
s q[11];
cx q[54], q[65];
t q[23];
t q[61];
h q[0];
cx q[45], q[32];
cx q[51], q[15];
cx q[59], q[32];
s q[6];
h q[16];
cx q[47], q[2];
cx q[62], q[47];
h q[49];
ccx q[56], q[53], q[27];
t q[28];
h q[35];
s q[43];
s q[54];
s q[13];
t q[23];
h q[6];
s q[34];
h q[34];
s q[65];
s q[66];
s q[55];
h q[50];
t q[7];
ccx q[30], q[13], q[12];
h q[37];
cx q[3], q[48];
ccx q[10], q[34], q[33];
h q[69];
ccx q[56], q[35], q[38];
h q[42];
cx q[2], q[17];
t q[39];
h q[63];
ccx q[67], q[18], q[42];
h q[30];
s q[53];
s q[1];
cx q[16], q[3];
h q[57];
t q[4];
s q[15];
cx q[39], q[19];
cx q[32], q[21];
cx q[25], q[30];
cx q[26], q[28];
ccx q[69], q[63], q[33];
t q[10];
s q[31];
ccx q[12], q[7], q[32];
ccx q[24], q[55], q[20];
h q[11];
t q[1];
cx q[0], q[1];
sdg q[0];
cx q[2], q[3];
tdg q[3];
sdg q[2];
cx q[2], q[3];
sdg q[5];
cx q[4], q[5];
sdg q[5];
cx q[4], q[5];
sdg q[4];
cx q[7], q[6];
tdg q[7];
tdg q[6];
sdg q[9];
cx q[8], q[9];
sdg q[8];
y q[10];
cx q[12], q[11];
sdg q[11];
tdg q[11];
cx q[14], q[13];
sdg q[13];
tdg q[14];
cx q[14], q[13];
sdg q[13];
cx q[15], q[16];
cx q[18], q[17];
tdg q[17];
cx q[19], q[20];
tdg q[20];
cx q[19], q[20];
tdg q[19];
sdg q[22];
sdg q[21];
cx q[21], q[22];
sdg q[21];
x q[23];
cx q[25], q[26];
sdg q[29];
sdg q[28];
cx q[30], q[31];
sdg q[32];
t q[54];
cx q[42], q[55];
s q[54];
h q[57];
t q[47];
cx q[66], q[58];
cx q[35], q[45];
t q[53];
h q[50];
h q[48];
h q[42];
s q[41];
h q[66];
ccx q[60], q[63], q[65];
s q[69];
cx q[55], q[59];
t q[35];
t q[66];
cx q[63], q[35];
cx q[53], q[41];
ccx q[39], q[35], q[46];
t q[57];
h q[37];
ccx q[48], q[38], q[41];
s q[53];
h q[57];
h q[64];
h q[37];
ccx q[63], q[47], q[69];
h q[38];
h q[45];
ccx q[45], q[54], q[36];
h q[52];
h q[37];
ccx q[50], q[37], q[51];
cx q[70], q[59];
cx q[71], q[44];
cx q[72], q[23];
cx q[73], q[16];
cx q[74], q[58];
cx q[75], q[62];
cx q[76], q[6];
cx q[77], q[64];
