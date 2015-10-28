@javac garch.java
@pause
@jar cvfm garch.jar mymanifest garch.class generalCarcass.class linearGarch.class nonLinearGarch.class commonVariables.class commonFunctions.class launcher.class gija.class
@echo -
@del *.class
@pause
