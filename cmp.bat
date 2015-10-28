@"C:\Program Files\Java\jdk1.6.0_18\bin\javac" garch.java
@pause
@"C:\Program Files\Java\jdk1.6.0_18\bin\jar" cvfm garch.jar mymanifest garch.class generalCarcass.class linearGarch.class nonLinearGarch.class commonVariables.class commonFunctions.class launcher.class gija.class
@echo -
@del *.class
@pause
