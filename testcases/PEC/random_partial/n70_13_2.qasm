OPENQASM 2.0;
include "qelib1.inc";
qreg q[81];
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
s q[11];
t q[25];
ccx q[65], q[18], q[9];
s q[9];
cx q[65], q[9];
cx q[44], q[66];
ccx q[45], q[28], q[42];
s q[5];
h q[65];
h q[55];
cx q[50], q[40];
t q[20];
t q[6];
t q[52];
ccx q[18], q[17], q[36];
t q[33];
t q[17];
t q[0];
t q[3];
s q[1];
ccx q[42], q[1], q[14];
t q[5];
ccx q[57], q[24], q[11];
s q[2];
h q[0];
cx q[35], q[0];
ccx q[57], q[17], q[18];
ccx q[59], q[24], q[6];
s q[69];
cx q[11], q[8];
t q[31];
t q[50];
h q[41];
t q[49];
t q[58];
t q[60];
h q[46];
ccx q[20], q[10], q[34];
t q[32];
ccx q[50], q[57], q[68];
s q[44];
s q[46];
cx q[21], q[41];
ccx q[65], q[35], q[56];
h q[12];
cx q[38], q[46];
t q[53];
ccx q[14], q[39], q[21];
t q[16];
t q[22];
h q[58];
s q[4];
t q[47];
ccx q[64], q[14], q[58];
cx q[32], q[40];
t q[41];
cx q[1], q[24];
t q[15];
cx q[19], q[26];
h q[4];
h q[34];
s q[67];
s q[13];
ccx q[41], q[26], q[11];
ccx q[43], q[33], q[29];
h q[40];
ccx q[65], q[9], q[59];
s q[67];
s q[0];
t q[45];
s q[62];
h q[16];
h q[57];
cx q[53], q[17];
cx q[4], q[35];
ccx q[56], q[69], q[60];
cx q[12], q[54];
cx q[47], q[33];
cx q[36], q[62];
h q[62];
ccx q[38], q[65], q[6];
s q[17];
t q[42];
h q[69];
t q[7];
h q[21];
h q[33];
s q[31];
h q[6];
t q[49];
cx q[2], q[67];
h q[8];
cx q[61], q[3];
cx q[31], q[65];
h q[53];
ccx q[38], q[37], q[40];
ccx q[65], q[62], q[43];
ccx q[60], q[36], q[32];
s q[51];
ccx q[20], q[59], q[36];
ccx q[17], q[32], q[6];
ccx q[50], q[1], q[52];
h q[35];
h q[1];
cx q[58], q[17];
cx q[38], q[32];
s q[5];
cx q[60], q[26];
cx q[16], q[36];
t q[17];
t q[1];
ccx q[12], q[0], q[35];
cx q[53], q[45];
cx q[36], q[8];
t q[66];
cx q[56], q[66];
h q[11];
s q[36];
h q[8];
s q[8];
cx q[10], q[48];
ccx q[58], q[18], q[9];
s q[4];
cx q[7], q[55];
h q[47];
t q[53];
t q[34];
h q[8];
ccx q[52], q[29], q[45];
ccx q[34], q[14], q[18];
ccx q[21], q[40], q[44];
h q[63];
s q[52];
t q[10];
ccx q[28], q[17], q[39];
cx q[57], q[17];
h q[48];
h q[31];
t q[63];
t q[13];
h q[12];
h q[34];
s q[62];
s q[69];
s q[28];
cx q[64], q[7];
s q[18];
cx q[4], q[46];
cx q[6], q[25];
h q[9];
t q[2];
t q[32];
h q[67];
s q[32];
ccx q[37], q[24], q[43];
ccx q[46], q[29], q[12];
cx q[68], q[61];
s q[40];
h q[39];
h q[50];
cx q[11], q[22];
h q[1];
s q[62];
ccx q[52], q[12], q[23];
s q[13];
cx q[64], q[6];
s q[66];
cx q[28], q[3];
cx q[18], q[22];
t q[0];
s q[45];
ccx q[15], q[60], q[30];
s q[2];
ccx q[19], q[22], q[23];
ccx q[15], q[48], q[20];
ccx q[39], q[4], q[58];
t q[46];
cx q[65], q[19];
h q[27];
s q[55];
h q[3];
t q[7];
cx q[30], q[32];
s q[24];
s q[11];
cx q[13], q[7];
t q[11];
ccx q[64], q[16], q[45];
t q[32];
t q[60];
h q[0];
ccx q[2], q[44], q[54];
ccx q[9], q[46], q[57];
ccx q[32], q[35], q[20];
t q[31];
cx q[61], q[24];
t q[46];
cx q[8], q[4];
t q[56];
h q[47];
h q[39];
s q[61];
s q[61];
cx q[3], q[65];
s q[58];
t q[65];
s q[12];
h q[63];
h q[27];
s q[21];
cx q[0], q[1];
sdg q[1];
cx q[0], q[1];
tdg q[0];
cx q[2], q[3];
tdg q[3];
tdg q[2];
cx q[2], q[3];
sdg q[2];
cx q[4], q[5];
tdg q[4];
sdg q[7];
cx q[10], q[9];
tdg q[9];
cx q[10], q[9];
tdg q[9];
tdg q[11];
cx q[11], q[12];
cx q[14], q[13];
tdg q[13];
cx q[14], q[13];
tdg q[13];
cx q[15], q[16];
tdg q[15];
cx q[18], q[17];
sdg q[17];
cx q[20], q[21];
sdg q[20];
cx q[22], q[23];
tdg q[22];
sdg q[22];
sdg q[24];
cx q[24], q[25];
tdg q[25];
cx q[24], q[25];
tdg q[24];
tdg q[27];
cx q[26], q[27];
tdg q[26];
cx q[28], q[29];
sdg q[28];
tdg q[29];
tdg q[28];
cx q[28], q[29];
sdg q[30];
sdg q[31];
cx q[30], q[31];
sdg q[30];
tdg q[32];
cx q[33], q[32];
sdg q[32];
y q[34];
ccx q[56], q[57], q[67];
h q[55];
s q[45];
ccx q[56], q[69], q[46];
s q[44];
t q[59];
t q[37];
cx q[59], q[69];
ccx q[61], q[55], q[39];
s q[62];
cx q[68], q[50];
h q[42];
s q[53];
ccx q[55], q[43], q[53];
ccx q[67], q[53], q[37];
cx q[49], q[35];
t q[57];
t q[58];
s q[66];
t q[40];
t q[42];
s q[57];
t q[46];
t q[54];
cx q[39], q[46];
cx q[57], q[66];
cx q[56], q[60];
s q[59];
h q[39];
ccx q[64], q[62], q[38];
t q[36];
h q[56];
cx q[55], q[40];
h q[55];
s q[48];
cx q[70], q[65];
cx q[71], q[39];
cx q[72], q[21];
cx q[73], q[60];
cx q[74], q[50];
cx q[75], q[13];
cx q[76], q[5];
cx q[77], q[26];
cx q[78], q[45];
cx q[79], q[16];
cx q[80], q[38];
