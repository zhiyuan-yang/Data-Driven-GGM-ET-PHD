for k=1:1:length(truth)
    hold on;
    plot(truth{k}.position(1,:),truth{k}.position(2,:),'bo');
    if ~isempty(output_GGIW_phd.state_X_hat{k})
        plot(output_GGIW_phd.state_X_hat{k}(1,:),output_GGIW_phd.state_X_hat{k}(2,:),'rx');
    end
    hold off;
    pause(0.25);
%     xlim([-500 500]);
%     ylim([-50 50]);
end