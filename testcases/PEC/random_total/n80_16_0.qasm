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
s q[75];
h q[69];
h q[72];
s q[62];
t q[8];
t q[9];
ccx q[53], q[13], q[65];
ccx q[13], q[59], q[22];
ccx q[21], q[18], q[34];
s q[32];
ccx q[49], q[35], q[20];
s q[27];
ccx q[26], q[3], q[47];
s q[13];
s q[28];
t q[11];
t q[74];
cx q[23], q[22];
s q[41];
ccx q[34], q[12], q[49];
ccx q[10], q[12], q[30];
cx q[12], q[5];
cx q[41], q[33];
t q[33];
ccx q[0], q[71], q[54];
h q[15];
h q[51];
t q[20];
ccx q[9], q[45], q[73];
s q[47];
h q[75];
s q[56];
cx q[78], q[22];
t q[7];
ccx q[33], q[65], q[49];
t q[61];
s q[39];
cx q[24], q[28];
h q[40];
h q[7];
t q[62];
cx q[55], q[26];
ccx q[44], q[28], q[40];
ccx q[67], q[13], q[36];
cx q[70], q[75];
ccx q[46], q[67], q[69];
ccx q[74], q[59], q[53];
ccx q[25], q[73], q[41];
ccx q[36], q[68], q[23];
h q[19];
ccx q[33], q[9], q[24];
ccx q[44], q[65], q[38];
cx q[79], q[42];
t q[21];
cx q[6], q[68];
cx q[16], q[79];
s q[78];
s q[40];
t q[25];
ccx q[73], q[60], q[1];
t q[78];
cx q[58], q[59];
s q[50];
s q[18];
t q[41];
cx q[40], q[72];
t q[56];
t q[31];
h q[11];
ccx q[28], q[58], q[56];
ccx q[37], q[44], q[29];
t q[28];
t q[7];
s q[56];
cx q[37], q[78];
t q[23];
h q[70];
ccx q[75], q[9], q[21];
s q[9];
s q[48];
h q[76];
t q[49];
t q[17];
s q[66];
h q[34];
cx q[24], q[30];
t q[30];
t q[39];
t q[4];
cx q[43], q[22];
t q[21];
t q[59];
s q[77];
ccx q[48], q[63], q[0];
h q[35];
ccx q[40], q[3], q[69];
cx q[69], q[61];
t q[69];
cx q[77], q[36];
cx q[30], q[16];
ccx q[65], q[53], q[61];
s q[50];
h q[42];
t q[76];
cx q[49], q[79];
t q[37];
t q[17];
t q[17];
s q[33];
s q[1];
t q[61];
h q[69];
h q[66];
ccx q[21], q[33], q[48];
h q[40];
h q[51];
ccx q[4], q[68], q[71];
h q[6];
h q[63];
ccx q[49], q[43], q[56];
s q[75];
t q[78];
t q[68];
t q[68];
s q[60];
ccx q[23], q[28], q[15];
ccx q[48], q[51], q[42];
t q[1];
ccx q[68], q[10], q[20];
h q[77];
h q[61];
h q[45];
h q[9];
cx q[52], q[4];
h q[74];
ccx q[47], q[54], q[59];
cx q[38], q[74];
ccx q[15], q[54], q[16];
t q[0];
t q[58];
h q[9];
cx q[71], q[9];
s q[5];
t q[37];
t q[50];
s q[50];
h q[26];
h q[9];
cx q[33], q[12];
ccx q[38], q[4], q[35];
cx q[45], q[0];
t q[24];
s q[36];
s q[26];
ccx q[42], q[43], q[23];
ccx q[1], q[53], q[2];
ccx q[56], q[8], q[1];
h q[8];
t q[22];
h q[63];
ccx q[18], q[2], q[74];
ccx q[46], q[27], q[18];
s q[52];
ccx q[41], q[29], q[79];
ccx q[64], q[46], q[58];
ccx q[74], q[50], q[10];
s q[18];
cx q[72], q[53];
h q[64];
t q[19];
s q[35];
ccx q[38], q[14], q[70];
t q[15];
t q[53];
ccx q[21], q[39], q[69];
cx q[51], q[70];
h q[12];
cx q[6], q[1];
s q[67];
t q[43];
t q[52];
t q[50];
t q[11];
s q[50];
ccx q[30], q[61], q[54];
s q[29];
t q[62];
s q[76];
ccx q[37], q[52], q[67];
h q[70];
t q[28];
cx q[9], q[66];
t q[55];
s q[8];
s q[58];
t q[37];
h q[72];
cx q[31], q[78];
s q[5];
ccx q[70], q[31], q[36];
ccx q[3], q[72], q[27];
cx q[34], q[48];
t q[66];
h q[7];
t q[54];
ccx q[22], q[43], q[10];
s q[61];
cx q[16], q[14];
t q[76];
t q[55];
cx q[55], q[64];
ccx q[30], q[24], q[68];
h q[40];
s q[69];
ccx q[0], q[38], q[24];
cx q[56], q[48];
cx q[76], q[26];
ccx q[64], q[33], q[30];
ccx q[13], q[0], q[36];
h q[29];
h q[76];
s q[63];
cx q[76], q[46];
ccx q[78], q[58], q[37];
ccx q[36], q[6], q[76];
s q[52];
h q[34];
h q[36];
t q[7];
ccx q[68], q[54], q[25];
cx q[46], q[76];
ccx q[49], q[42], q[21];
t q[16];
cx q[64], q[47];
h q[18];
t q[47];
ccx q[34], q[38], q[58];
cx q[56], q[79];
ccx q[7], q[44], q[75];
ccx q[46], q[29], q[31];
