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
ccx q[27], q[29], q[17];
ccx q[3], q[38], q[21];
cx q[38], q[21];
s q[22];
t q[34];
cx q[30], q[13];
cx q[23], q[38];
cx q[25], q[12];
s q[25];
s q[15];
s q[19];
h q[3];
cx q[35], q[29];
h q[21];
h q[18];
ccx q[13], q[12], q[14];
h q[17];
ccx q[1], q[4], q[3];
s q[19];
cx q[6], q[27];
s q[29];
h q[11];
ccx q[19], q[39], q[37];
s q[16];
s q[1];
ccx q[23], q[28], q[1];
h q[3];
ccx q[21], q[6], q[32];
t q[31];
cx q[24], q[16];
t q[24];
s q[36];
h q[31];
s q[25];
t q[19];
ccx q[17], q[39], q[8];
cx q[32], q[4];
ccx q[3], q[23], q[30];
t q[14];
h q[28];
ccx q[10], q[38], q[34];
s q[10];
t q[21];
h q[31];
cx q[12], q[4];
s q[37];
s q[26];
s q[31];
h q[25];
s q[13];
ccx q[33], q[4], q[25];
h q[35];
ccx q[1], q[10], q[32];
cx q[19], q[18];
ccx q[35], q[30], q[7];
ccx q[1], q[12], q[38];
t q[3];
s q[33];
ccx q[28], q[38], q[4];
ccx q[20], q[32], q[22];
t q[17];
ccx q[23], q[28], q[35];
s q[27];
h q[21];
t q[26];
h q[12];
s q[30];
ccx q[28], q[19], q[35];
t q[34];
s q[2];
s q[18];
t q[5];
cx q[8], q[27];
ccx q[5], q[12], q[4];
ccx q[11], q[4], q[17];
cx q[13], q[32];
ccx q[5], q[20], q[7];
h q[37];
s q[9];
h q[8];
s q[11];
t q[11];
s q[36];
ccx q[2], q[38], q[11];
t q[20];
cx q[32], q[37];
ccx q[11], q[36], q[20];
ccx q[18], q[9], q[25];
cx q[7], q[11];
t q[10];
t q[28];
t q[15];
h q[28];
cx q[0], q[19];
s q[32];
h q[38];
s q[22];
s q[11];
t q[5];
t q[10];
ccx q[37], q[5], q[6];
h q[35];
h q[36];
s q[25];
s q[35];
t q[6];
h q[39];
ccx q[21], q[34], q[38];
cx q[37], q[24];
s q[29];
h q[4];
h q[1];
cx q[38], q[22];
h q[18];
h q[4];
s q[12];
s q[34];
h q[31];
h q[12];
s q[29];