OPENQASM 2.0;
include "qelib1.inc";
qreg q[80];
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
h q[14];
t q[5];
s q[61];
ccx q[43], q[18], q[68];
h q[10];
s q[0];
h q[15];
cx q[3], q[9];
ccx q[42], q[45], q[23];
cx q[47], q[12];
h q[38];
s q[57];
t q[9];
cx q[21], q[65];
h q[8];
cx q[32], q[30];
s q[8];
cx q[62], q[61];
cx q[50], q[3];
s q[16];
h q[15];
s q[50];
ccx q[50], q[63], q[14];
cx q[43], q[25];
s q[48];
t q[58];
ccx q[67], q[34], q[33];
h q[19];
t q[41];
cx q[6], q[0];
s q[42];
s q[38];
s q[20];
t q[68];
t q[4];
h q[14];
s q[4];
cx q[31], q[40];
ccx q[39], q[38], q[16];
t q[13];
h q[69];
ccx q[17], q[38], q[59];
s q[15];
ccx q[66], q[11], q[54];
h q[20];
s q[31];
cx q[47], q[41];
s q[33];
cx q[47], q[23];
s q[44];
ccx q[48], q[36], q[7];
t q[14];
cx q[45], q[33];
t q[11];
t q[17];
ccx q[56], q[20], q[37];
h q[28];
cx q[63], q[33];
cx q[59], q[2];
ccx q[5], q[39], q[52];
cx q[67], q[6];
s q[56];
h q[8];
s q[16];
cx q[9], q[3];
h q[0];
ccx q[47], q[36], q[16];
t q[12];
h q[53];
h q[59];
cx q[55], q[60];
s q[14];
ccx q[3], q[61], q[48];
t q[32];
t q[56];
t q[39];
h q[27];
cx q[11], q[46];
s q[58];
ccx q[2], q[43], q[13];
t q[4];
h q[20];
t q[51];
h q[46];
s q[22];
s q[19];
ccx q[40], q[67], q[7];
h q[29];
ccx q[29], q[47], q[17];
cx q[2], q[54];
h q[15];
ccx q[50], q[53], q[23];
s q[31];
t q[18];
s q[31];
s q[67];
cx q[32], q[30];
s q[15];
s q[62];
t q[43];
h q[11];
t q[5];
s q[1];
h q[69];
ccx q[47], q[58], q[57];
ccx q[19], q[15], q[32];
ccx q[10], q[52], q[69];
ccx q[59], q[37], q[55];
ccx q[6], q[68], q[56];
h q[41];
h q[48];
s q[43];
h q[16];
h q[69];
t q[68];
s q[1];
cx q[8], q[19];
h q[45];
cx q[21], q[25];
h q[55];
ccx q[23], q[27], q[12];
t q[4];
ccx q[35], q[27], q[55];
s q[68];
s q[36];
h q[9];
ccx q[59], q[52], q[49];
s q[43];
h q[9];
h q[24];
s q[20];
t q[0];
s q[56];
s q[28];
h q[55];
ccx q[47], q[34], q[53];
t q[7];
ccx q[40], q[20], q[24];
cx q[1], q[41];
ccx q[66], q[28], q[36];
s q[16];
cx q[54], q[23];
ccx q[20], q[63], q[23];
ccx q[11], q[39], q[49];
cx q[49], q[42];
s q[44];
t q[36];
cx q[49], q[29];
cx q[60], q[68];
t q[41];
ccx q[23], q[36], q[42];
s q[31];
t q[61];
s q[38];
ccx q[13], q[22], q[33];
cx q[23], q[18];
ccx q[3], q[9], q[4];
t q[53];
h q[61];
ccx q[64], q[47], q[11];
cx q[58], q[25];
ccx q[68], q[1], q[57];
s q[21];
ccx q[42], q[48], q[19];
h q[57];
t q[66];
h q[63];
s q[1];
ccx q[27], q[25], q[39];
cx q[36], q[50];
cx q[15], q[32];
s q[62];
ccx q[59], q[64], q[38];
h q[8];
s q[31];
h q[43];
cx q[30], q[16];
cx q[41], q[45];
cx q[15], q[24];
t q[28];
t q[34];
h q[20];
h q[67];
h q[20];
cx q[4], q[23];
t q[58];
h q[0];
h q[54];
cx q[25], q[40];
t q[46];
s q[17];
s q[40];
t q[55];
s q[25];
s q[60];
h q[69];
ccx q[47], q[40], q[28];
cx q[18], q[7];
t q[50];
cx q[14], q[56];
s q[48];
h q[53];
ccx q[60], q[59], q[15];
t q[55];
h q[47];
h q[60];
t q[11];
s q[2];
cx q[35], q[45];
t q[23];
x q[0];
z q[2];
sdg q[3];
cx q[5], q[6];
tdg q[5];
tdg q[6];
sdg q[5];
cx q[5], q[6];
cx q[7], q[8];
tdg q[8];
sdg q[7];
cx q[7], q[8];
tdg q[9];
sdg q[9];
cx q[9], q[10];
cx q[12], q[11];
tdg q[11];
cx q[12], q[11];
x q[11];
tdg q[14];
cx q[13], q[14];
tdg q[13];
x q[15];
x q[16];
tdg q[18];
cx q[17], q[18];
sdg q[18];
cx q[17], q[18];
sdg q[17];
cx q[20], q[19];
x q[19];
cx q[22], q[21];
sdg q[21];
cx q[24], q[23];
sdg q[23];
z q[25];
cx q[26], q[27];
tdg q[27];
sdg q[26];
cx q[26], q[27];
tdg q[26];
x q[28];
sdg q[31];
cx q[32], q[31];
x q[31];
x q[33];
x q[34];
ccx q[57], q[56], q[62];
ccx q[49], q[58], q[51];
h q[68];
ccx q[68], q[55], q[49];
s q[61];
t q[59];
t q[67];
cx q[41], q[57];
cx q[52], q[36];
t q[39];
h q[49];
s q[62];
cx q[35], q[36];
t q[37];
s q[43];
h q[63];
s q[44];
s q[45];
cx q[59], q[54];
ccx q[68], q[52], q[39];
s q[53];
ccx q[55], q[65], q[41];
ccx q[52], q[51], q[50];
t q[49];
ccx q[55], q[39], q[61];
s q[62];
h q[69];
s q[42];
s q[49];
s q[60];
cx q[56], q[57];
t q[64];
t q[57];
t q[40];
ccx q[65], q[42], q[53];
cx q[70], q[49];
cx q[71], q[63];
cx q[72], q[32];
cx q[73], q[48];
cx q[74], q[65];
cx q[75], q[31];
cx q[76], q[19];
cx q[77], q[62];
cx q[78], q[17];
cx q[79], q[34];
