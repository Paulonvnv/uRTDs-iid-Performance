hmean<-function(x) {1/mean(1/x)}
sd.hmean<-function(x) {sqrt((mean(1/x))^(-4)*var(1/x)/length(x))}
