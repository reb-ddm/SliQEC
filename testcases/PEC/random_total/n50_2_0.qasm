OPENQASM 2.0;
include "qelib1.inc";
qreg q[50];
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
s q[26];
t q[23];
ccx q[43], q[38], q[21];
h q[48];
s q[47];
t q[9];
t q[10];
t q[18];
cx q[32], q[1];
h q[21];
cx q[42], q[37];
h q[39];
s q[21];
s q[20];
ccx q[35], q[17], q[4];
ccx q[5], q[9], q[20];
ccx q[10], q[22], q[36];
s q[1];
t q[39];
t q[20];
s q[5];
t q[6];
t q[2];
h q[33];
cx q[1], q[23];
ccx q[37], q[6], q[49];
t q[25];
t q[37];
ccx q[44], q[10], q[6];
s q[44];
t q[17];
h q[41];
t q[45];
t q[31];
ccx q[43], q[21], q[32];
s q[14];
cx q[12], q[23];
ccx q[1], q[4], q[7];
h q[24];
ccx q[14], q[40], q[31];
h q[25];
cx q[14], q[0];
ccx q[25], q[11], q[33];
ccx q[23], q[48], q[21];
s q[45];
t q[2];
ccx q[41], q[5], q[20];
cx q[2], q[10];
t q[4];
cx q[49], q[30];
s q[37];
t q[24];
s q[34];
h q[11];
t q[44];
s q[48];
cx q[0], q[15];
cx q[37], q[5];
h q[36];
s q[37];
s q[39];
ccx q[34], q[23], q[10];
s q[1];
cx q[35], q[39];
cx q[34], q[23];
ccx q[22], q[42], q[30];
h q[27];
cx q[28], q[38];
cx q[25], q[39];
s q[25];
cx q[2], q[24];
t q[2];
t q[13];
h q[38];
t q[19];
t q[39];
t q[12];
s q[9];
t q[36];
ccx q[26], q[4], q[25];
s q[45];
cx q[16], q[45];
h q[35];
t q[0];
cx q[3], q[0];
t q[23];
h q[16];
s q[48];
t q[49];
h q[17];
t q[20];
cx q[28], q[9];
t q[14];
s q[41];
ccx q[33], q[19], q[42];
h q[49];
s q[17];
cx q[16], q[46];
h q[21];
s q[29];
h q[40];
cx q[7], q[35];
t q[36];
s q[9];
h q[21];
h q[17];
ccx q[48], q[2], q[44];
h q[0];
s q[12];
h q[8];
h q[39];
h q[36];
ccx q[30], q[0], q[34];
h q[26];
cx q[28], q[5];
cx q[14], q[28];
h q[40];
ccx q[29], q[20], q[19];
cx q[32], q[3];
ccx q[23], q[39], q[30];
cx q[14], q[39];
s q[46];
h q[19];
t q[2];
t q[5];
t q[31];
ccx q[10], q[40], q[14];
ccx q[2], q[25], q[43];
s q[16];
t q[32];
s q[14];
ccx q[34], q[5], q[12];
ccx q[13], q[47], q[39];
s q[8];
s q[8];
ccx q[13], q[24], q[33];
cx q[33], q[19];
s q[3];
cx q[8], q[22];
h q[26];
h q[34];
ccx q[4], q[36], q[40];
ccx q[20], q[45], q[15];
ccx q[48], q[14], q[43];
h q[2];
h q[43];
t q[26];
t q[6];
h q[1];
h q[3];
