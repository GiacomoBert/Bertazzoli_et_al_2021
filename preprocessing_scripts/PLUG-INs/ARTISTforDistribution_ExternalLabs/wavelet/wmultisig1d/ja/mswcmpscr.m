% MSWCMPSCR  1 ���������M���̃E�F�[�u���b�g���k�X�R�A
%
%   [THR,L2SCR,N0SCR,IDXSORT] = MSWCMPSCR(DEC) �́A4 �̍s����v�Z���܂��B
%   �X���b�V���z�[���h THR�A���k�X�R�A L2SCR �� N0SCR �ƃC���f�b�N�X IDXSORT�B
%   ���� DEC �́Adetail �� (�I�v�V�����Ƃ���) approximation �W���̘A����
%   ������E�F�[�u���b�g�W�� CFS �̍s��ɑΉ����܂��B
%       CFS = [cd{DEC.LEV} , ... , cd{1}] �܂���
%       CFS = [ca , cd{DEC.LEV} , ... , cd{1}]
%   �A���́ADEC.dirDec �� 'r' (�܂��� 'c') �̏ꍇ�A�s���� (�܂��͗����) 
%   �ɍs�Ȃ��܂��B
%
%   �I���W�i���M���̐� NbSIG �Ɗe�M���ɑ΂���W���̐� NbCFS �� (���ׂāA
%   �܂��� detail �W���̂�) �^����ƁACFS �� NbSIG �s NbCFS ��̍s���
%   �Ȃ�܂��B
%	������:
%     - THR, L2SCR, N0SCR �́ANbSIG �s (NbCFS+1) ��̍s��ł��B
%     - IDXSORT �� NbSIG �s NbCFS ��̍s��ł��B
%   ����:
%     - THR(:,2:end) �́A��Βl�ɑ΂��ď����̍s�Ŋi�[���ꂽ CFS �Ɠ�����
%       �Ȃ�܂��B
%       �܂��ATHR(:,1) = 0 �ł��B
%     i �Ԗڂ̐M���ɑ΂��āA
%     - L2SCR(i,j) �́ACFS(i,j-1) �ɓ������X���b�V���z�[���h�ɑΉ�����
%       �ۑ����ꂽ�G�l���M�[ (L2-�m����) �̊����ł� 
%       (2 <= j <= NbCFS)�B�܂��AL2SCR(:,1) = 100 �ł��B
%     - N0SCR(i,j) �́ACFS(i,j-1) �ɓ������X���b�V���z�[���h�ɑΉ�����
%       0 �̊����ł� (2 <= j <= NbCFS)�B�܂��AN0SCR(:,1) = 0 �ł��B
%
%   ����� 3 �̃I�v�V�������͂��g�p�ł��܂��B
%     [...] = MSWCMPSCR(...,S_or_H,KEEPAPP,IDXSIG)
%       - S_or_H  ('s' �܂��� 'h') �́A�\�t�g�A�܂��̓n�[�h�X���b�V��
%         �z�[���h���Ƃ����Ӗ��ł� (�ڍׂɂ��Ă� MSWTHRESH ���Q��)�B
%       - KEEPAPP (true �܂��� false)�BKEEPAPP �� true �ɓ������ꍇ�A
%         approximation �W���͕ێ�����܂��B
%       - IDXSIG �́A�����M���̃C���f�b�N�X���܂ރx�N�g���A�܂���
%         ������ 'all' �ł��B
%   �f�t�H���g�́A���ꂼ����̂悤�ɂȂ�܂�: 'h', false, 'all'.
%
%   �Q�l mdwtdec, mdwtrec, ddencmp, wdencmp


%   Copyright 1995-2007 The MathWorks, Inc.
