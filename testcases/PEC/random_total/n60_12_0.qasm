OPENQASM 2.0;
include "qelib1.inc";
qreg q[60];
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
t q[12];
ccx q[31], q[0], q[28];
h q[17];
cx q[12], q[6];
ccx q[49], q[26], q[34];
t q[11];
h q[59];
t q[7];
ccx q[7], q[11], q[9];
t q[31];
cx q[51], q[10];
cx q[18], q[42];
s q[1];
cx q[18], q[9];
s q[48];
h q[53];
ccx q[41], q[20], q[47];
h q[58];
h q[44];
s q[33];
h q[52];
s q[15];
s q[24];
t q[21];
h q[26];
t q[37];
ccx q[23], q[33], q[55];
s q[19];
cx q[43], q[50];
ccx q[13], q[57], q[22];
s q[26];
s q[21];
ccx q[44], q[50], q[29];
h q[12];
s q[5];
ccx q[22], q[19], q[48];
s q[12];
h q[31];
h q[48];
ccx q[6], q[44], q[19];
ccx q[9], q[39], q[30];
ccx q[42], q[28], q[26];
h q[0];
ccx q[59], q[50], q[54];
h q[35];
s q[53];
s q[58];
t q[44];
ccx q[32], q[15], q[26];
s q[25];
ccx q[54], q[45], q[46];
t q[54];
h q[44];
s q[40];
h q[15];
ccx q[38], q[21], q[1];
cx q[58], q[2];
s q[4];
ccx q[41], q[31], q[20];
cx q[9], q[59];
s q[48];
h q[17];
ccx q[14], q[52], q[45];
cx q[35], q[51];
ccx q[39], q[43], q[49];
s q[19];
h q[27];
cx q[22], q[45];
ccx q[58], q[50], q[2];
s q[50];
s q[51];
t q[9];
ccx q[51], q[12], q[39];
h q[43];
h q[3];
h q[41];
s q[49];
h q[52];
s q[13];
s q[41];
ccx q[2], q[38], q[59];
ccx q[45], q[16], q[41];
h q[34];
t q[12];
t q[26];
s q[44];
t q[9];
h q[35];
cx q[15], q[56];
h q[17];
cx q[11], q[8];
ccx q[30], q[50], q[7];
ccx q[6], q[58], q[3];
ccx q[32], q[21], q[20];
ccx q[19], q[57], q[6];
s q[45];
h q[33];
t q[25];
cx q[41], q[10];
cx q[50], q[40];
s q[50];
cx q[11], q[14];
ccx q[12], q[56], q[59];
t q[36];
t q[15];
ccx q[53], q[0], q[45];
h q[24];
t q[58];
t q[25];
s q[51];
s q[37];
s q[42];
cx q[50], q[13];
ccx q[5], q[16], q[3];
cx q[9], q[28];
h q[59];
cx q[27], q[10];
t q[23];
t q[40];
t q[44];
cx q[32], q[29];
t q[51];
ccx q[12], q[5], q[46];
ccx q[11], q[32], q[47];
ccx q[45], q[16], q[42];
s q[1];
cx q[5], q[39];
t q[36];
ccx q[31], q[52], q[19];
cx q[25], q[5];
t q[21];
cx q[19], q[55];
t q[34];
h q[2];
t q[9];
cx q[32], q[27];
h q[43];
s q[28];
s q[59];
cx q[21], q[51];
cx q[19], q[3];
t q[11];
ccx q[16], q[37], q[7];
h q[14];
h q[28];
t q[12];
s q[40];
s q[51];
s q[12];
cx q[2], q[16];
ccx q[28], q[19], q[15];
s q[41];
h q[10];
t q[24];
ccx q[30], q[18], q[21];
t q[17];
s q[52];
h q[9];
t q[37];
h q[27];
s q[47];
t q[56];
h q[45];
cx q[10], q[44];
h q[42];
cx q[57], q[48];
ccx q[7], q[55], q[22];
s q[22];
ccx q[55], q[39], q[32];
h q[36];
t q[40];
t q[7];
ccx q[39], q[5], q[41];
cx q[34], q[32];
h q[31];
s q[8];
ccx q[23], q[5], q[10];
t q[33];
s q[45];
h q[35];
