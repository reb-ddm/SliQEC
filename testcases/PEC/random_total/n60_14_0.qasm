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
cx q[46], q[9];
cx q[23], q[48];
t q[18];
ccx q[43], q[48], q[22];
h q[11];
s q[56];
ccx q[55], q[5], q[23];
cx q[54], q[47];
h q[53];
t q[8];
s q[10];
s q[36];
ccx q[46], q[0], q[17];
h q[2];
ccx q[11], q[29], q[51];
t q[46];
t q[25];
t q[57];
h q[42];
s q[12];
h q[10];
s q[53];
s q[2];
s q[39];
ccx q[22], q[37], q[49];
ccx q[20], q[12], q[51];
h q[6];
h q[8];
h q[54];
ccx q[55], q[43], q[6];
h q[27];
cx q[12], q[34];
ccx q[7], q[9], q[59];
ccx q[54], q[18], q[20];
ccx q[37], q[14], q[55];
s q[27];
h q[30];
s q[7];
cx q[47], q[17];
ccx q[16], q[22], q[27];
t q[9];
ccx q[53], q[57], q[39];
cx q[10], q[50];
h q[35];
h q[34];
cx q[57], q[22];
cx q[50], q[19];
h q[1];
s q[11];
h q[26];
cx q[1], q[5];
ccx q[43], q[56], q[51];
h q[5];
ccx q[8], q[42], q[43];
s q[20];
ccx q[17], q[19], q[25];
t q[34];
ccx q[50], q[20], q[22];
s q[20];
t q[53];
h q[48];
t q[18];
h q[51];
cx q[57], q[36];
ccx q[20], q[17], q[31];
s q[45];
s q[55];
ccx q[1], q[39], q[11];
h q[58];
t q[34];
t q[1];
s q[20];
h q[56];
t q[58];
cx q[48], q[39];
ccx q[58], q[40], q[22];
ccx q[39], q[31], q[20];
t q[12];
h q[36];
cx q[16], q[5];
t q[32];
h q[1];
cx q[29], q[12];
t q[47];
cx q[46], q[9];
t q[29];
h q[15];
h q[33];
h q[38];
s q[38];
ccx q[50], q[17], q[8];
cx q[44], q[54];
t q[0];
s q[55];
cx q[38], q[54];
cx q[29], q[4];
ccx q[57], q[35], q[13];
h q[9];
t q[43];
t q[29];
cx q[26], q[20];
h q[53];
cx q[50], q[37];
h q[25];
s q[9];
h q[55];
h q[37];
t q[46];
h q[16];
ccx q[13], q[14], q[58];
t q[24];
h q[56];
h q[52];
ccx q[18], q[35], q[52];
s q[41];
ccx q[35], q[44], q[37];
s q[35];
s q[33];
ccx q[42], q[55], q[51];
h q[15];
h q[38];
h q[18];
ccx q[11], q[59], q[42];
h q[45];
h q[5];
s q[15];
t q[49];
cx q[44], q[25];
t q[59];
s q[2];
t q[17];
t q[26];
cx q[19], q[57];
t q[27];
h q[56];
cx q[35], q[50];
ccx q[54], q[15], q[17];
h q[10];
s q[50];
h q[36];
ccx q[16], q[24], q[35];
t q[46];
h q[38];
ccx q[43], q[55], q[16];
cx q[18], q[16];
s q[37];
h q[51];
t q[52];
h q[13];
ccx q[27], q[5], q[28];
t q[40];
h q[8];
h q[3];
t q[27];
cx q[43], q[26];
ccx q[51], q[2], q[56];
s q[29];
ccx q[4], q[51], q[7];
ccx q[57], q[49], q[44];
ccx q[45], q[24], q[17];
ccx q[5], q[34], q[58];
ccx q[47], q[13], q[1];
t q[18];
ccx q[29], q[39], q[9];
ccx q[22], q[58], q[55];
cx q[13], q[40];
s q[59];
cx q[54], q[42];
ccx q[48], q[7], q[38];
t q[56];
t q[14];
t q[12];
ccx q[12], q[55], q[13];
cx q[44], q[28];
ccx q[30], q[16], q[41];
s q[37];
ccx q[56], q[22], q[58];
cx q[2], q[57];
s q[2];
ccx q[32], q[13], q[24];
