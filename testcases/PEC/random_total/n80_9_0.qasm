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
t q[33];
h q[41];
cx q[34], q[48];
s q[0];
t q[19];
h q[71];
cx q[52], q[61];
ccx q[64], q[20], q[36];
t q[64];
s q[44];
cx q[41], q[69];
cx q[1], q[7];
cx q[33], q[46];
s q[73];
h q[31];
cx q[47], q[18];
cx q[62], q[60];
cx q[73], q[5];
h q[67];
t q[14];
h q[9];
h q[72];
ccx q[27], q[61], q[12];
h q[29];
cx q[35], q[40];
ccx q[54], q[2], q[50];
s q[44];
t q[16];
t q[75];
h q[78];
s q[56];
s q[26];
s q[47];
h q[30];
s q[25];
cx q[16], q[57];
ccx q[75], q[49], q[7];
h q[30];
s q[11];
t q[50];
t q[32];
h q[9];
t q[58];
t q[7];
h q[47];
cx q[30], q[11];
cx q[10], q[47];
t q[26];
ccx q[16], q[52], q[54];
cx q[20], q[56];
cx q[33], q[1];
t q[51];
cx q[51], q[0];
t q[65];
cx q[11], q[10];
s q[27];
t q[10];
ccx q[10], q[46], q[42];
h q[47];
cx q[53], q[22];
s q[63];
h q[38];
s q[16];
cx q[65], q[47];
cx q[36], q[12];
t q[68];
cx q[23], q[78];
t q[68];
ccx q[41], q[78], q[30];
s q[41];
s q[37];
h q[38];
t q[31];
cx q[8], q[39];
ccx q[41], q[38], q[53];
ccx q[43], q[9], q[12];
s q[44];
cx q[23], q[66];
t q[37];
t q[47];
s q[66];
h q[36];
s q[1];
cx q[10], q[43];
s q[0];
t q[10];
s q[41];
ccx q[59], q[42], q[24];
cx q[7], q[25];
ccx q[57], q[69], q[51];
s q[68];
cx q[65], q[20];
cx q[48], q[12];
cx q[62], q[69];
h q[39];
h q[52];
ccx q[30], q[57], q[19];
s q[73];
h q[13];
cx q[31], q[76];
s q[53];
ccx q[54], q[71], q[65];
t q[21];
s q[40];
t q[56];
h q[70];
t q[57];
t q[0];
ccx q[71], q[16], q[1];
t q[29];
ccx q[44], q[68], q[57];
s q[2];
cx q[12], q[23];
ccx q[48], q[59], q[38];
t q[72];
h q[7];
cx q[60], q[71];
t q[29];
s q[43];
ccx q[7], q[11], q[2];
cx q[29], q[6];
s q[62];
h q[48];
cx q[15], q[12];
h q[5];
t q[14];
s q[51];
h q[35];
t q[65];
cx q[57], q[36];
s q[42];
h q[70];
cx q[25], q[54];
t q[24];
ccx q[39], q[60], q[59];
cx q[4], q[18];
h q[3];
t q[68];
cx q[1], q[51];
h q[19];
t q[71];
t q[69];
s q[51];
h q[22];
s q[47];
h q[43];
t q[8];
t q[26];
cx q[6], q[27];
ccx q[63], q[32], q[9];
s q[72];
ccx q[21], q[75], q[13];
h q[27];
cx q[66], q[24];
ccx q[54], q[33], q[41];
s q[28];
s q[5];
t q[72];
cx q[35], q[49];
cx q[58], q[35];
ccx q[44], q[13], q[71];
s q[18];
ccx q[25], q[69], q[49];
t q[75];
h q[75];
h q[54];
cx q[76], q[69];
t q[62];
t q[33];
t q[19];
ccx q[12], q[78], q[41];
ccx q[65], q[44], q[8];
s q[0];
h q[35];
s q[62];
cx q[36], q[39];
h q[23];
cx q[17], q[56];
s q[58];
ccx q[15], q[28], q[19];
ccx q[75], q[25], q[17];
h q[6];
s q[41];
ccx q[52], q[66], q[44];
cx q[60], q[55];
cx q[11], q[9];
h q[59];
cx q[14], q[39];
t q[65];
cx q[12], q[8];
cx q[68], q[64];
ccx q[69], q[39], q[59];
t q[13];
h q[28];
s q[57];
h q[63];
cx q[28], q[9];
h q[5];
h q[9];
cx q[56], q[10];
t q[5];
cx q[50], q[46];
s q[49];
h q[21];
t q[26];
h q[50];
cx q[37], q[1];
s q[21];
t q[49];
ccx q[5], q[22], q[0];
t q[27];
s q[76];
s q[24];
s q[31];
h q[67];
ccx q[36], q[27], q[39];
h q[17];
ccx q[3], q[4], q[53];
cx q[25], q[67];
h q[26];
ccx q[35], q[2], q[77];
ccx q[9], q[78], q[54];
h q[4];
cx q[50], q[66];
cx q[45], q[78];
h q[7];
t q[6];
cx q[76], q[7];
ccx q[62], q[76], q[74];
cx q[64], q[60];
s q[38];
h q[20];
h q[5];
h q[62];
cx q[21], q[1];
cx q[52], q[70];
cx q[54], q[22];
h q[60];
cx q[0], q[22];
h q[27];
