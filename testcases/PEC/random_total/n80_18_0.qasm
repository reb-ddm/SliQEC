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
t q[57];
t q[56];
cx q[65], q[68];
h q[13];
t q[43];
h q[66];
ccx q[47], q[68], q[9];
s q[11];
ccx q[23], q[2], q[10];
cx q[68], q[18];
ccx q[33], q[28], q[75];
cx q[49], q[9];
ccx q[6], q[46], q[47];
ccx q[68], q[9], q[61];
t q[16];
ccx q[19], q[27], q[73];
h q[45];
ccx q[50], q[37], q[53];
ccx q[64], q[13], q[0];
cx q[57], q[10];
h q[29];
ccx q[37], q[46], q[58];
ccx q[71], q[31], q[16];
cx q[55], q[69];
ccx q[27], q[74], q[28];
ccx q[8], q[58], q[21];
s q[30];
cx q[38], q[41];
ccx q[37], q[10], q[24];
ccx q[54], q[20], q[16];
h q[15];
t q[11];
ccx q[15], q[78], q[77];
t q[22];
h q[50];
s q[30];
t q[66];
h q[7];
h q[46];
ccx q[78], q[15], q[1];
cx q[45], q[3];
cx q[73], q[1];
h q[28];
cx q[48], q[30];
t q[38];
s q[19];
s q[47];
cx q[9], q[0];
t q[56];
cx q[5], q[55];
cx q[48], q[50];
h q[59];
h q[59];
s q[43];
h q[31];
ccx q[58], q[62], q[11];
t q[76];
t q[34];
cx q[10], q[6];
ccx q[27], q[9], q[56];
ccx q[21], q[9], q[32];
t q[76];
h q[75];
cx q[29], q[10];
ccx q[47], q[37], q[69];
cx q[61], q[44];
cx q[15], q[57];
t q[70];
h q[68];
h q[25];
t q[63];
h q[79];
ccx q[36], q[15], q[26];
s q[13];
h q[50];
t q[63];
ccx q[70], q[7], q[78];
t q[29];
h q[17];
cx q[29], q[58];
cx q[79], q[7];
s q[26];
t q[73];
t q[39];
s q[8];
ccx q[21], q[23], q[58];
ccx q[5], q[77], q[56];
t q[27];
h q[6];
t q[74];
h q[11];
cx q[41], q[22];
s q[29];
ccx q[7], q[57], q[3];
t q[76];
s q[19];
cx q[39], q[16];
t q[33];
h q[21];
h q[60];
t q[28];
h q[7];
h q[6];
cx q[21], q[11];
cx q[70], q[9];
t q[11];
t q[19];
t q[18];
ccx q[67], q[19], q[79];
h q[5];
h q[11];
t q[46];
cx q[72], q[50];
ccx q[48], q[41], q[13];
t q[8];
h q[29];
s q[16];
s q[38];
cx q[30], q[11];
s q[32];
t q[35];
t q[35];
t q[0];
ccx q[12], q[14], q[35];
h q[53];
s q[0];
cx q[22], q[66];
cx q[25], q[47];
h q[71];
t q[71];
s q[53];
t q[1];
h q[13];
ccx q[13], q[32], q[6];
ccx q[15], q[19], q[45];
cx q[18], q[26];
t q[25];
s q[78];
s q[28];
t q[71];
s q[24];
h q[17];
t q[46];
h q[69];
t q[45];
t q[79];
s q[49];
t q[78];
s q[9];
cx q[46], q[25];
cx q[67], q[45];
ccx q[39], q[3], q[58];
t q[55];
h q[22];
h q[15];
h q[74];
s q[19];
h q[37];
ccx q[75], q[46], q[58];
cx q[46], q[62];
ccx q[4], q[29], q[0];
s q[16];
cx q[52], q[43];
h q[40];
ccx q[22], q[12], q[38];
cx q[64], q[2];
t q[1];
h q[41];
s q[56];
s q[42];
cx q[9], q[2];
s q[71];
ccx q[36], q[44], q[78];
h q[68];
t q[57];
ccx q[15], q[79], q[7];
s q[19];
h q[43];
s q[55];
cx q[69], q[77];
cx q[0], q[51];
ccx q[61], q[46], q[64];
h q[33];
s q[53];
s q[1];
cx q[63], q[7];
h q[31];
cx q[70], q[79];
cx q[31], q[46];
s q[21];
t q[33];
cx q[50], q[35];
t q[58];
s q[16];
t q[59];
t q[60];
t q[72];
cx q[62], q[51];
ccx q[13], q[14], q[49];
s q[6];
s q[61];
t q[51];
h q[62];
cx q[76], q[66];
t q[60];
s q[28];
t q[30];
h q[44];
ccx q[33], q[12], q[34];
cx q[45], q[7];
cx q[54], q[20];
h q[61];
h q[21];
cx q[50], q[62];
s q[74];
t q[46];
h q[47];
h q[47];
ccx q[69], q[20], q[1];
s q[31];
t q[14];
cx q[27], q[44];
cx q[10], q[50];
h q[60];
h q[8];
h q[41];
s q[40];
s q[55];
ccx q[26], q[50], q[10];
t q[56];
s q[72];
h q[45];
ccx q[25], q[37], q[28];
s q[55];
h q[60];
h q[47];
cx q[27], q[41];
h q[8];
s q[0];
cx q[2], q[66];