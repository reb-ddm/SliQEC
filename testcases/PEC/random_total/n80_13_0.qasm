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
t q[19];
s q[34];
s q[18];
cx q[16], q[78];
cx q[40], q[70];
h q[2];
s q[73];
h q[58];
cx q[56], q[65];
cx q[38], q[2];
ccx q[7], q[2], q[77];
s q[29];
cx q[3], q[13];
h q[54];
t q[65];
ccx q[18], q[11], q[46];
cx q[35], q[4];
h q[16];
cx q[17], q[67];
h q[25];
cx q[47], q[45];
h q[40];
s q[34];
h q[1];
t q[3];
ccx q[71], q[60], q[66];
s q[3];
h q[3];
t q[50];
ccx q[40], q[22], q[68];
h q[45];
t q[10];
cx q[64], q[26];
t q[39];
cx q[32], q[31];
cx q[27], q[47];
cx q[48], q[78];
s q[61];
ccx q[23], q[1], q[74];
h q[62];
s q[59];
ccx q[14], q[64], q[39];
s q[34];
ccx q[63], q[38], q[3];
t q[15];
t q[0];
t q[6];
s q[24];
t q[65];
t q[36];
ccx q[79], q[15], q[18];
ccx q[39], q[8], q[73];
t q[9];
t q[27];
t q[8];
cx q[36], q[65];
t q[44];
cx q[49], q[62];
cx q[53], q[26];
cx q[71], q[20];
t q[73];
ccx q[28], q[71], q[66];
cx q[4], q[59];
s q[68];
h q[13];
t q[47];
t q[73];
h q[55];
ccx q[0], q[5], q[68];
cx q[76], q[21];
cx q[61], q[68];
ccx q[14], q[50], q[12];
h q[54];
ccx q[14], q[72], q[70];
ccx q[36], q[7], q[79];
s q[50];
h q[30];
h q[38];
ccx q[62], q[66], q[30];
s q[68];
s q[31];
ccx q[40], q[73], q[75];
t q[23];
s q[5];
h q[48];
t q[73];
h q[19];
s q[7];
cx q[35], q[68];
t q[61];
h q[66];
t q[79];
t q[57];
h q[78];
h q[33];
t q[24];
s q[20];
ccx q[33], q[46], q[23];
cx q[75], q[7];
t q[38];
s q[67];
h q[9];
cx q[77], q[5];
h q[27];
cx q[15], q[62];
s q[15];
s q[68];
h q[35];
cx q[57], q[28];
ccx q[65], q[57], q[61];
h q[53];
s q[49];
s q[26];
s q[7];
s q[18];
h q[35];
cx q[20], q[75];
t q[6];
t q[49];
ccx q[35], q[46], q[67];
ccx q[19], q[28], q[40];
t q[73];
ccx q[44], q[5], q[54];
t q[75];
s q[74];
s q[49];
cx q[23], q[34];
t q[54];
s q[74];
cx q[5], q[42];
s q[55];
h q[71];
t q[69];
ccx q[61], q[24], q[58];
t q[61];
t q[66];
cx q[66], q[9];
t q[7];
s q[10];
s q[59];
cx q[14], q[59];
s q[55];
s q[56];
h q[59];
s q[42];
cx q[2], q[42];
cx q[10], q[67];
t q[14];
s q[46];
cx q[27], q[49];
cx q[56], q[58];
ccx q[16], q[10], q[40];
ccx q[48], q[21], q[38];
ccx q[21], q[24], q[65];
ccx q[54], q[73], q[14];
s q[45];
ccx q[70], q[75], q[30];
ccx q[4], q[27], q[66];
t q[6];
s q[23];
ccx q[47], q[9], q[7];
h q[43];
h q[77];
cx q[77], q[31];
ccx q[59], q[26], q[14];
s q[22];
t q[68];
t q[53];
ccx q[63], q[15], q[78];
t q[57];
cx q[49], q[27];
h q[17];
t q[58];
cx q[21], q[24];
s q[53];
cx q[69], q[40];
t q[34];
ccx q[76], q[4], q[42];
ccx q[58], q[8], q[77];
h q[57];
t q[35];
cx q[33], q[29];
h q[59];
h q[12];
ccx q[40], q[69], q[27];
h q[44];
s q[78];
s q[18];
t q[8];
t q[48];
t q[53];
t q[29];
t q[5];
h q[7];
h q[64];
cx q[7], q[28];
cx q[64], q[1];
ccx q[15], q[59], q[17];
t q[61];
ccx q[58], q[62], q[35];
ccx q[15], q[33], q[31];
ccx q[66], q[39], q[22];
cx q[65], q[45];
ccx q[11], q[57], q[68];
cx q[66], q[32];
t q[74];
h q[0];
h q[33];
t q[33];
cx q[35], q[69];
h q[4];
s q[27];
ccx q[57], q[25], q[16];
h q[4];
cx q[17], q[15];
s q[16];
h q[11];
cx q[52], q[70];
s q[12];
ccx q[42], q[57], q[24];
t q[39];
ccx q[41], q[61], q[16];
h q[45];
s q[41];
ccx q[47], q[23], q[63];
ccx q[77], q[70], q[65];
h q[60];
s q[48];
ccx q[28], q[66], q[6];
t q[9];
s q[51];
h q[15];
ccx q[20], q[21], q[63];
s q[69];
s q[57];
ccx q[16], q[77], q[18];
h q[73];
h q[75];
ccx q[55], q[57], q[1];
ccx q[22], q[63], q[25];