OPENQASM 2.0;
include "qelib1.inc";
qreg q[40];
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
t q[0];
ccx q[9], q[16], q[14];
t q[19];
s q[15];
t q[9];
ccx q[5], q[6], q[31];
t q[39];
h q[14];
t q[10];
h q[17];
h q[1];
h q[19];
t q[20];
cx q[26], q[8];
ccx q[33], q[7], q[38];
cx q[28], q[22];
t q[36];
cx q[6], q[28];
cx q[7], q[27];
h q[19];
cx q[16], q[10];
cx q[38], q[8];
ccx q[26], q[12], q[8];
s q[33];
h q[14];
s q[23];
t q[18];
cx q[31], q[32];
ccx q[2], q[21], q[17];
t q[27];
ccx q[31], q[37], q[1];
t q[27];
h q[12];
h q[15];
t q[32];
s q[1];
cx q[5], q[37];
ccx q[14], q[5], q[16];
h q[37];
cx q[24], q[7];
cx q[5], q[0];
t q[10];
t q[31];
s q[13];
t q[22];
s q[9];
s q[31];
s q[31];
s q[9];
t q[17];
h q[9];
t q[4];
cx q[32], q[34];
s q[27];
s q[28];
ccx q[3], q[23], q[6];
cx q[17], q[36];
ccx q[4], q[31], q[22];
s q[5];
s q[4];
h q[27];
h q[21];
ccx q[20], q[7], q[10];
t q[19];
ccx q[7], q[34], q[8];
t q[26];
s q[9];
h q[21];
t q[12];
s q[18];
cx q[0], q[33];
t q[34];
h q[39];
h q[13];
s q[36];
cx q[1], q[8];
cx q[25], q[28];
ccx q[29], q[20], q[33];
cx q[30], q[5];
t q[4];
cx q[12], q[5];
h q[22];
h q[24];
cx q[26], q[12];
t q[31];
h q[20];
ccx q[8], q[19], q[18];
ccx q[0], q[8], q[3];
s q[35];
t q[20];
cx q[38], q[16];
cx q[14], q[12];
t q[30];
t q[11];
h q[17];
s q[17];
s q[32];
ccx q[25], q[7], q[26];
cx q[28], q[26];
cx q[5], q[0];
s q[19];
h q[14];
ccx q[8], q[11], q[24];
s q[12];
cx q[39], q[33];
h q[11];
h q[2];
ccx q[8], q[17], q[38];
cx q[12], q[4];
h q[1];
t q[13];
ccx q[33], q[1], q[20];
ccx q[4], q[34], q[20];
s q[39];
ccx q[4], q[28], q[30];
cx q[35], q[36];
ccx q[28], q[36], q[29];
s q[37];
h q[34];
ccx q[3], q[38], q[2];
