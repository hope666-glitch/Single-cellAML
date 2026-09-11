pdf(
  file = "AUCAML94.pdf",  
  width = 4.5,  
  height = 4.5,
  paper = "special"
)
plot.roc(roc1, col = "#6699CC", lty=1,lwd=6,legacy.axes = TRUE,
         main="",
         xlab = "1-Specificity", ylab = "Sensitivity",
         xlim = c(1, 0), ylim = c(0, 1),
         font.lab = 2,
         cex.lab = 1.2, 
         cex.axis = 1.2,
         font = 2)
plot.roc(roc2, col = "#D40D8C", add = T,lty=1,lwd=6)
plot.roc(roc3, col = "#FBB01A", add = T,lty=1,lwd=6)

text(0.05,0.25,paste("AUC = ",round(roc1$auc,3)),
     cex = 1.2, font = 2,
     pos = 2)
text(0.05,0.15,paste("AUC = ",round(roc2$auc,3)),
     cex = 1.2, font = 2,
     pos = 2)
text(0.05,0.05,paste("AUC = ",round(roc3$auc,3)),
     cex = 1.2, font = 2,
     pos = 2)
legend(x=0.12,y=0.35,c("","",""),col=c("#0070B8","#D40D8C","#FBB01A"),
       lty = 1,
       lwd = 5,
       y.intersp = 1,    
       cex = 1.2,         
       text.font = 2, 
       bty = "n",          
       bg = "transparent")
dev.off()
plot.roc(roc2, col = "#D40D8C", lty=1,lwd=6,legacy.axes = TRUE,
         main="",
         xlab = "1-Specificity", ylab = "Sensitivity",
         xlim = c(1, 0), ylim = c(0, 1),
         font.lab = 2,
         cex.lab = 1.2, 
         cex.axis = 1.2,
         font = 2)
text(0.05,0.25,paste("AUC = ",round(roc2$auc,3)),
     cex = 1.2, font = 2,
     pos = 2)
legend(x=0.12,y=0.35,c(""),col=c("#D40D8C"),
       lty = 1,
       lwd = 5,
       y.intersp = 1,    
       cex = 1.2,         
       text.font = 2, 
       bty = "n",          
       bg = "transparent")
plot.roc(roc2, col = "#FBB01A", lty=1,lwd=6,legacy.axes = TRUE,
         main="",
         xlab = "1-Specificity", ylab = "Sensitivity",
         xlim = c(1, 0), ylim = c(0, 1),
         font.lab = 2,
         cex.lab = 1.2, 
         cex.axis = 1.2,
         font = 2)
text(0.05,0.25,paste("AUC = ",round(roc2$auc,3)),
     cex = 1.2, font = 2,
     pos = 2)
legend(x=0.12,y=0.35,c(""),col=c("#FBB01A"),
       lty = 1,
       lwd = 5,
       y.intersp = 1,    
       cex = 1.2,         
       text.font = 2, 
       bty = "n",          
       bg = "transparent")
dev.off()
full_df$risk_eln <- predict(fit2, type = "response")
full_df$risk_combo <- predict(fit3, type = "response")

dca_res <- dca(fit2, fit3, model.names = c("Model1","Model2"))
dca_res <- dca(fit2, fit3,
               model.names = c("Model1","Model2"),
               threshold_range = c(0,0.75))
p_dca <- ggplot(dca_res, linetype = F)+
  geom_line(aes(color = model), linewidth = 1.5)+
  theme_bw(base_size = 12)+
  labs(x = "Threshold probability", y = "Net benefit")+
  theme(legend.position = "bottomright")
p_dca
ggplot(dca_res, linetype = F) +  
  geom_line(aes(color = model), linewidth = 1.5) +  
  
  scale_color_manual(
    values = c(
      "Model2" = "#E53E3E",       
      "Model1" = "#2196F3",  
      "All" = "#805AD5",   
      "None" = "#f9a826"             
    )
  ) +
  
  guides(linetype = "none") +
  
  labs(
    title = "",
    x = "Threshold Probability (%)",
    y = "Net Benefit"
  ) +
  
  theme_minimal() +
  theme(
    
    plot.title = element_text(
      size = 18,          
      face = "bold",      
      hjust = 1,          
      vjust = 1,          
      margin = margin(b = -20) 
    ),
    
    axis.title.x = element_text(
      size = 16,          
      face = "bold",      
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 16,          
      face = "bold",      
      margin = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      size = 16,          
      face = "bold",      
      color = "black"     
    ),
    axis.text.y = element_text(
      size = 16,          
      face = "bold",      
      color = "black"     
    ),
    
    legend.text = element_text(
      size = 13,          
      face = "bold"       
    ),
    legend.title = element_blank(), 
    legend.position = "right",      
    legend.key.size = unit(0.6, "cm"),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black",    
                             linewidth = 1.2),
    axis.ticks = element_line(
      color = "black",    
      linewidth = 1.2     
    ),
    
    axis.ticks.length = unit(0.25, "cm"), 
    
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )
ggsave(
  "dca_style_optimized.pdf",
  width = 6, height =4 , dpi = 300, device = "pdf"
)