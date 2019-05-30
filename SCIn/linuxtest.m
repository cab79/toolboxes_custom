lj = Labjack_linux('verbose',true) %open labJack verbosely
%lj.toggleFIO(9) %toggle FIO1 between low and high
%lj.timedTTL(9,1) % send a 200ms timed TTL pulse on FIO3
%lj.strobeWord
lj.setDIOValue(5)
pause(0.1)
lj.setDIOValue(0)
%out=lj.rawWrite(00000001)