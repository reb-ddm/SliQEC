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
h q[22];
h q[2];
h q[35];
s q[19];
s q[30];
s q[13];
s q[20];
ccx q[23], q[4], q[16];
ccx q[32], q[31], q[7];
h q[7];
cx q[15], q[19];
ccx q[11], q[17], q[30];
cx q[13], q[31];
ccx q[28], q[24], q[0];
s q[23];
ccx q[39], q[17], q[6];
s q[1];
t q[11];
ccx q[36], q[19], q[11];
s q[22];
t q[2];
h q[29];
ccx q[3], q[0], q[5];
ccx q[0], q[14], q[8];
ccx q[1], q[22], q[10];
h q[7];
ccx q[23], q[7], q[35];
s q[33];
ccx q[21], q[14], q[2];
h q[19];
cx q[30], q[3];
h q[27];
ccx q[23], q[26], q[11];
h q[12];
ccx q[13], q[31], q[11];
ccx q[7], q[6], q[28];
s q[1];
ccx q[9], q[12], q[16];
t q[14];
cx q[37], q[20];
t q[9];
h q[24];
h q[33];
cx q[14], q[3];
h q[25];
h q[36];
ccx q[31], q[32], q[21];
ccx q[29], q[23], q[20];
s q[2];
h q[13];
t q[20];
s q[35];
cx q[21], q[10];
s q[37];
s q[39];
t q[13];
cx q[22], q[27];
t q[23];
t q[0];
h q[2];
cx q[37], q[32];
cx q[32], q[35];
s q[15];
s q[18];
s q[7];
ccx q[0], q[30], q[29];
t q[25];
t q[26];
cx q[27], q[22];
h q[33];
h q[16];
t q[5];
ccx q[27], q[39], q[29];
t q[34];
ccx q[34], q[18], q[35];
cx q[36], q[15];
h q[30];
h q[19];
cx q[10], q[20];
t q[21];
ccx q[19], q[38], q[12];
s q[37];
cx q[9], q[14];
h q[0];
h q[39];
h q[5];
s q[31];
t q[2];
ccx q[8], q[26], q[5];
ccx q[19], q[16], q[39];
ccx q[30], q[12], q[5];
ccx q[7], q[24], q[6];
t q[9];
ccx q[12], q[26], q[15];
s q[36];
h q[29];
ccx q[13], q[39], q[35];
h q[38];
cx q[37], q[26];
h q[23];
ccx q[4], q[25], q[1];
s q[8];
h q[1];
cx q[31], q[26];
ccx q[20], q[37], q[9];
s q[5];
ccx q[24], q[36], q[30];
ccx q[21], q[24], q[7];
ccx q[39], q[19], q[23];
h q[10];
h q[1];
t q[36];
ccx q[20], q[7], q[36];
ccx q[15], q[6], q[29];
cx q[31], q[14];
h q[34];
s q[18];
t q[14];
t q[25];
h q[12];