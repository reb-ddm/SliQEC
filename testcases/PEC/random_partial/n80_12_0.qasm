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
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
h q[36];
s q[75];
s q[39];
s q[59];
h q[44];
s q[75];
ccx q[10], q[52], q[55];
t q[58];
h q[64];
cx q[17], q[73];
cx q[25], q[10];
ccx q[21], q[63], q[53];
t q[66];
t q[62];
s q[6];
ccx q[6], q[36], q[20];
t q[50];
ccx q[77], q[65], q[24];
h q[22];
cx q[52], q[18];
s q[65];
s q[73];
s q[1];
ccx q[17], q[14], q[40];
s q[57];
t q[53];
ccx q[52], q[24], q[59];
t q[68];
cx q[53], q[66];
s q[6];
cx q[7], q[10];
s q[33];
cx q[37], q[34];
cx q[67], q[31];
h q[20];
cx q[68], q[46];
h q[49];
t q[11];
ccx q[19], q[70], q[17];
s q[0];
cx q[74], q[13];
s q[0];
ccx q[78], q[12], q[57];
ccx q[42], q[32], q[77];
h q[7];
h q[15];
ccx q[2], q[5], q[43];
s q[57];
t q[52];
s q[43];
h q[16];
cx q[56], q[78];
h q[65];
cx q[9], q[57];
s q[41];
cx q[17], q[70];
t q[12];
ccx q[64], q[77], q[78];
h q[62];
cx q[38], q[28];
t q[60];
s q[40];
t q[72];
ccx q[10], q[22], q[23];
h q[49];
s q[45];
h q[21];
h q[57];
ccx q[64], q[3], q[6];
s q[75];
cx q[62], q[74];
s q[40];
h q[35];
h q[47];
ccx q[61], q[18], q[72];
h q[59];
s q[24];
cx q[63], q[49];
t q[39];
t q[13];
cx q[14], q[39];
t q[52];
cx q[76], q[54];
t q[14];
t q[1];
t q[3];
s q[62];
ccx q[25], q[63], q[62];
cx q[42], q[8];
s q[29];
h q[15];
t q[55];
s q[73];
t q[69];
s q[26];
cx q[11], q[22];
s q[63];
cx q[1], q[7];
ccx q[78], q[39], q[46];
cx q[57], q[18];
s q[77];
t q[64];
t q[7];
cx q[22], q[56];
ccx q[56], q[55], q[16];
s q[9];
cx q[34], q[13];
ccx q[61], q[74], q[53];
s q[9];
s q[7];
t q[15];
ccx q[34], q[33], q[60];
ccx q[21], q[17], q[35];
s q[45];
h q[6];
ccx q[27], q[29], q[49];
t q[22];
t q[28];
h q[0];
t q[14];
s q[69];
t q[13];
ccx q[38], q[33], q[71];
t q[57];
t q[28];
s q[2];
s q[45];
cx q[66], q[54];
ccx q[35], q[28], q[17];
h q[48];
cx q[32], q[20];
ccx q[4], q[40], q[37];
s q[69];
cx q[50], q[74];
s q[56];
h q[73];
h q[42];
t q[33];
h q[45];
ccx q[50], q[51], q[59];
s q[12];
h q[67];
cx q[73], q[75];
ccx q[9], q[51], q[38];
h q[63];
t q[49];
ccx q[13], q[51], q[78];
ccx q[33], q[35], q[57];
t q[26];
h q[55];
h q[12];
s q[59];
ccx q[47], q[76], q[30];
s q[58];
ccx q[69], q[53], q[32];
h q[25];
cx q[10], q[47];
h q[26];
ccx q[60], q[31], q[69];
t q[54];
s q[51];
h q[79];
h q[60];
cx q[10], q[49];
cx q[13], q[16];
ccx q[69], q[3], q[61];
ccx q[15], q[14], q[42];
ccx q[24], q[72], q[79];
cx q[68], q[22];
s q[32];
t q[64];
s q[28];
h q[11];
h q[16];
t q[30];
cx q[5], q[66];
t q[39];
ccx q[23], q[60], q[74];
h q[28];
cx q[56], q[36];
ccx q[75], q[43], q[70];
cx q[73], q[79];
h q[23];
h q[70];
t q[15];
s q[42];
ccx q[48], q[1], q[57];
s q[22];
s q[30];
h q[76];
cx q[40], q[5];
s q[40];
s q[56];
cx q[17], q[8];
ccx q[71], q[30], q[36];
cx q[3], q[26];
s q[48];
h q[55];
ccx q[53], q[41], q[4];
t q[32];
s q[27];
s q[75];
ccx q[14], q[37], q[69];
s q[54];
h q[52];
t q[8];
h q[10];
cx q[42], q[63];
h q[24];
ccx q[2], q[17], q[64];
ccx q[68], q[69], q[11];
t q[38];
t q[38];
s q[29];
cx q[10], q[25];
t q[11];
t q[6];
s q[44];
t q[55];
h q[20];
t q[50];
h q[33];
h q[31];
s q[64];
ccx q[40], q[7], q[63];
s q[37];
cx q[70], q[52];
s q[2];
t q[15];
t q[25];
t q[16];
ccx q[8], q[44], q[59];
t q[66];
h q[36];
t q[4];
ccx q[15], q[28], q[47];
h q[59];
h q[56];
t q[5];
s q[33];
x q[0];
y q[1];
cx q[2], q[3];
tdg q[3];
tdg q[2];
cx q[2], q[3];
tdg q[5];
cx q[5], q[6];
sdg q[5];
cx q[8], q[7];
sdg q[7];
cx q[9], q[10];
sdg q[10];
tdg q[9];
sdg q[11];
cx q[12], q[11];
tdg q[11];
sdg q[13];
cx q[13], q[14];
tdg q[13];
x q[15];
tdg q[17];
sdg q[16];
cx q[16], q[17];
cx q[19], q[18];
x q[18];
cx q[20], q[21];
tdg q[20];
tdg q[22];
sdg q[25];
sdg q[24];
cx q[24], q[25];
sdg q[27];
sdg q[26];
cx q[27], q[26];
tdg q[26];
sdg q[28];
tdg q[30];
sdg q[32];
cx q[33], q[32];
tdg q[32];
cx q[35], q[34];
tdg q[34];
sdg q[37];
cx q[78], q[45];
ccx q[68], q[57], q[69];
cx q[69], q[78];
h q[65];
cx q[51], q[45];
ccx q[40], q[76], q[58];
h q[53];
h q[62];
cx q[68], q[79];
h q[43];
s q[46];
s q[65];
s q[64];
s q[57];
h q[43];
h q[70];
t q[55];
h q[72];
cx q[77], q[42];
ccx q[76], q[78], q[69];
s q[73];
ccx q[67], q[77], q[69];
t q[51];
ccx q[58], q[63], q[67];
ccx q[61], q[54], q[44];
h q[76];
s q[78];
ccx q[48], q[69], q[70];
cx q[61], q[58];
t q[77];
s q[78];
h q[78];
h q[77];
s q[74];
t q[65];
t q[59];
t q[50];
h q[77];
t q[77];
t q[52];
