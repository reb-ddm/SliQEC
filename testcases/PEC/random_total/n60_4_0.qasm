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
ccx q[40], q[52], q[12];
t q[16];
ccx q[16], q[35], q[28];
ccx q[33], q[23], q[44];
t q[47];
s q[34];
h q[36];
ccx q[43], q[36], q[7];
cx q[20], q[17];
cx q[38], q[30];
cx q[20], q[1];
t q[31];
t q[23];
ccx q[45], q[53], q[58];
cx q[53], q[21];
cx q[33], q[10];
t q[5];
t q[55];
s q[47];
t q[49];
ccx q[7], q[28], q[41];
t q[42];
s q[18];
h q[53];
ccx q[4], q[18], q[53];
cx q[57], q[48];
h q[47];
cx q[57], q[37];
cx q[11], q[3];
t q[50];
cx q[50], q[14];
cx q[37], q[58];
cx q[2], q[14];
t q[3];
s q[29];
h q[2];
ccx q[29], q[12], q[22];
t q[43];
ccx q[16], q[50], q[7];
ccx q[7], q[56], q[46];
h q[39];
ccx q[58], q[39], q[23];
h q[34];
cx q[55], q[56];
ccx q[47], q[43], q[29];
ccx q[0], q[44], q[17];
ccx q[54], q[13], q[56];
s q[48];
h q[47];
cx q[4], q[14];
cx q[39], q[44];
t q[28];
cx q[36], q[55];
h q[23];
h q[23];
s q[33];
ccx q[3], q[24], q[25];
cx q[36], q[33];
s q[42];
t q[50];
ccx q[59], q[12], q[33];
cx q[48], q[49];
cx q[26], q[45];
cx q[15], q[41];
s q[27];
h q[25];
t q[37];
ccx q[11], q[13], q[49];
t q[6];
t q[26];
h q[41];
t q[2];
ccx q[43], q[12], q[26];
h q[33];
cx q[53], q[25];
t q[55];
ccx q[50], q[47], q[2];
ccx q[12], q[35], q[57];
t q[41];
s q[29];
cx q[19], q[50];
s q[19];
ccx q[40], q[0], q[43];
s q[13];
cx q[19], q[33];
s q[50];
cx q[2], q[31];
s q[46];
h q[13];
cx q[59], q[7];
t q[31];
ccx q[33], q[21], q[58];
ccx q[46], q[41], q[18];
t q[29];
s q[44];
t q[0];
s q[43];
t q[52];
ccx q[23], q[28], q[48];
ccx q[39], q[36], q[42];
s q[21];
cx q[33], q[19];
ccx q[58], q[24], q[35];
cx q[17], q[16];
t q[9];
ccx q[54], q[5], q[3];
t q[20];
t q[18];
ccx q[45], q[31], q[12];
s q[29];
t q[1];
s q[1];
t q[51];
cx q[3], q[29];
t q[7];
cx q[10], q[45];
ccx q[19], q[43], q[35];
h q[20];
ccx q[58], q[30], q[6];
h q[33];
cx q[56], q[22];
s q[32];
cx q[57], q[2];
ccx q[22], q[23], q[55];
ccx q[6], q[11], q[0];
h q[0];
t q[58];
t q[59];
ccx q[40], q[41], q[47];
s q[33];
s q[31];
cx q[45], q[33];
t q[39];
h q[32];
cx q[45], q[48];
t q[54];
ccx q[56], q[26], q[46];
cx q[1], q[31];
h q[32];
ccx q[43], q[29], q[57];
ccx q[11], q[4], q[26];
t q[24];
h q[45];
s q[55];
s q[40];
ccx q[13], q[52], q[16];
ccx q[59], q[11], q[17];
h q[23];
ccx q[40], q[17], q[32];
t q[32];
t q[44];
cx q[54], q[39];
ccx q[16], q[6], q[58];
h q[48];
h q[33];
ccx q[19], q[11], q[16];
ccx q[54], q[59], q[8];
t q[21];
ccx q[11], q[30], q[48];
h q[9];
ccx q[14], q[38], q[15];
h q[33];
cx q[24], q[23];
cx q[47], q[43];
s q[37];
cx q[44], q[47];
ccx q[46], q[26], q[54];
h q[39];
s q[39];
h q[23];
t q[47];
t q[57];
cx q[37], q[1];
ccx q[4], q[8], q[47];
cx q[33], q[12];
t q[58];
cx q[16], q[51];
ccx q[46], q[20], q[51];
cx q[0], q[25];