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
t q[18];
ccx q[65], q[25], q[21];
h q[64];
s q[63];
ccx q[72], q[18], q[65];
s q[52];
ccx q[14], q[40], q[74];
cx q[71], q[54];
h q[54];
t q[38];
ccx q[65], q[14], q[38];
t q[5];
h q[54];
ccx q[32], q[51], q[0];
t q[71];
cx q[59], q[44];
cx q[59], q[28];
s q[32];
cx q[1], q[59];
ccx q[40], q[12], q[61];
s q[3];
t q[22];
h q[28];
ccx q[49], q[40], q[17];
t q[14];
cx q[38], q[54];
cx q[18], q[15];
t q[38];
h q[62];
cx q[70], q[74];
s q[17];
cx q[35], q[58];
s q[34];
t q[15];
h q[34];
ccx q[25], q[24], q[32];
s q[16];
h q[59];
ccx q[75], q[49], q[44];
cx q[49], q[8];
t q[59];
cx q[2], q[39];
t q[2];
ccx q[59], q[2], q[36];
cx q[52], q[43];
h q[69];
t q[43];
s q[50];
ccx q[19], q[25], q[21];
cx q[75], q[24];
ccx q[46], q[7], q[48];
h q[62];
s q[38];
h q[7];
h q[72];
h q[71];
ccx q[46], q[42], q[19];
s q[30];
s q[36];
t q[61];
cx q[14], q[1];
ccx q[27], q[25], q[19];
h q[49];
t q[58];
cx q[26], q[29];
cx q[7], q[75];
t q[75];
h q[13];
h q[20];
t q[65];
ccx q[61], q[32], q[59];
ccx q[40], q[47], q[6];
s q[8];
t q[57];
s q[41];
ccx q[59], q[40], q[37];
h q[27];
s q[64];
cx q[57], q[13];
cx q[44], q[23];
s q[7];
t q[18];
h q[15];
t q[8];
cx q[75], q[54];
ccx q[32], q[69], q[11];
t q[64];
cx q[39], q[10];
t q[62];
t q[71];
cx q[20], q[39];
ccx q[66], q[55], q[31];
cx q[79], q[49];
s q[21];
h q[11];
s q[24];
t q[79];
h q[65];
ccx q[64], q[5], q[59];
h q[51];
t q[42];
ccx q[35], q[44], q[62];
t q[73];
ccx q[79], q[54], q[12];
h q[37];
s q[54];
ccx q[60], q[11], q[26];
h q[73];
cx q[59], q[3];
t q[6];
h q[61];
ccx q[2], q[57], q[63];
s q[1];
cx q[42], q[2];
t q[15];
s q[31];
t q[60];
t q[9];
cx q[47], q[33];
t q[38];
t q[67];
ccx q[6], q[0], q[22];
t q[61];
t q[59];
t q[49];
t q[41];
s q[27];
h q[51];
h q[74];
s q[25];
s q[65];
h q[69];
t q[61];
cx q[66], q[16];
s q[21];
cx q[60], q[73];
h q[21];
ccx q[34], q[66], q[74];
h q[57];
t q[24];
t q[50];
t q[38];
s q[24];
h q[26];
ccx q[22], q[45], q[50];
t q[74];
cx q[18], q[66];
h q[56];
h q[43];
s q[77];
cx q[52], q[61];
h q[27];
cx q[34], q[11];
cx q[67], q[62];
cx q[51], q[27];
t q[8];
s q[48];
ccx q[7], q[32], q[49];
t q[59];
t q[44];
cx q[74], q[15];
s q[30];
ccx q[66], q[56], q[59];
s q[72];
t q[4];
s q[52];
ccx q[64], q[79], q[72];
h q[55];
s q[56];
h q[49];
cx q[43], q[5];
h q[8];
h q[65];
s q[32];
ccx q[19], q[21], q[14];
t q[45];
t q[48];
t q[26];
cx q[55], q[35];
h q[30];
ccx q[13], q[62], q[49];
cx q[67], q[9];
ccx q[44], q[1], q[19];
cx q[74], q[59];
cx q[55], q[51];
s q[11];
t q[28];
cx q[11], q[35];
ccx q[28], q[65], q[18];
cx q[39], q[9];
t q[73];
cx q[13], q[33];
ccx q[46], q[65], q[66];
ccx q[0], q[3], q[41];
ccx q[59], q[43], q[78];
cx q[7], q[59];
s q[51];
s q[29];
s q[4];
s q[76];
t q[58];
t q[61];
s q[22];
h q[27];
ccx q[56], q[12], q[77];
t q[21];
ccx q[51], q[8], q[7];
h q[49];
s q[22];
h q[73];
h q[11];
t q[18];
cx q[10], q[6];
h q[15];
t q[58];
cx q[62], q[2];
ccx q[45], q[77], q[19];
t q[3];
ccx q[49], q[3], q[39];
t q[38];
s q[65];
t q[54];
t q[75];
s q[77];
h q[21];
ccx q[24], q[27], q[13];
cx q[23], q[72];
h q[44];
t q[52];
h q[66];
ccx q[53], q[26], q[58];
ccx q[78], q[3], q[67];
s q[7];
cx q[10], q[13];
h q[6];
cx q[5], q[66];
ccx q[76], q[49], q[71];
cx q[69], q[58];
t q[49];
ccx q[49], q[60], q[13];
cx q[1], q[0];
sdg q[0];
cx q[2], q[3];
tdg q[3];
sdg q[2];
cx q[4], q[5];
tdg q[5];
sdg q[4];
cx q[4], q[5];
sdg q[4];
x q[6];
z q[8];
cx q[11], q[10];
tdg q[10];
cx q[11], q[10];
tdg q[10];
sdg q[13];
sdg q[16];
cx q[15], q[16];
tdg q[15];
sdg q[15];
cx q[18], q[17];
sdg q[17];
cx q[18], q[17];
sdg q[17];
z q[20];
tdg q[21];
cx q[22], q[21];
tdg q[21];
sdg q[21];
cx q[23], q[24];
sdg q[24];
sdg q[23];
tdg q[25];
cx q[26], q[25];
tdg q[25];
sdg q[28];
tdg q[27];
tdg q[29];
cx q[30], q[29];
sdg q[29];
cx q[31], q[32];
sdg q[31];
cx q[34], q[33];
sdg q[34];
tdg q[33];
cx q[34], q[33];
sdg q[33];
cx q[36], q[35];
tdg q[35];
tdg q[37];
s q[51];
s q[63];
h q[44];
cx q[74], q[72];
s q[46];
cx q[66], q[53];
h q[47];
s q[42];
cx q[46], q[73];
t q[65];
s q[56];
h q[63];
h q[51];
t q[75];
s q[47];
cx q[64], q[54];
cx q[58], q[67];
t q[71];
ccx q[48], q[51], q[42];
cx q[54], q[50];
h q[42];
ccx q[64], q[47], q[79];
ccx q[61], q[73], q[77];
s q[57];
t q[71];
ccx q[78], q[41], q[43];
cx q[43], q[66];
t q[73];
s q[67];
s q[57];
t q[51];
t q[50];
s q[77];
h q[52];
ccx q[69], q[76], q[58];
h q[71];
h q[61];
t q[70];
t q[55];
h q[75];
