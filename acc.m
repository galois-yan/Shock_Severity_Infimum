classdef acc < timeseries
    %
    %ACC A collection of methods dealing with acceleration time history
    %signal of mechanical shocks. This class is written as a subclass of
    %'timeseries' class, so all 'timeseries' methods can also be used.
    %
    %ACC Properties:
    %   Sf - Sample rate
    %   Time - Time column
    %   Data - Measured acceleration data column
    %   Length - Length of time series
    %
    %ACC Methods:
    %   acc - Constructor method to creat a ACC object.
    %   resample1 - Resample a time series.
    %   srs - Plot the SRS curve of the time series.
    %   bandpass - A bandpass filter.
    %   plot - An overload plot function for ACC object.
    %   fft - An overload fast fourier transform function for ACC object.
    %   fit - Shock waveform decomposition method.
    %   cwt - An overload continues wavelet transform plot.
    %   dwt - An overload discrete wavelet transform plot.
    %   sound - Play the time series as a sound, at specific sample rate.
    %   audiowrite - Write the time series as a sound in .wav format.
    %   cumtraapz - Overload numerical integration.
    %   diff - Overload numerical difference.
    %   extend - Extending the time series for a certain period.
    %   subplot - Plot all time series in a subplot view.
    %
    properties
        Sf  % Sample rate
    end

    methods
        function Acc = acc(t, y)
            %Acc = ACC(Sign) construct an ACC object from the measurement
            %matrix 'Sign'. The first column of 'Sign' shall be time
            %column, and the rest columns shall be acceleration data.

            Acc@timeseries(y, t-t(1));
            Acc.Sf = 1/(t(2)-t(1));
            %Acc=setuniformtime(Acc,'StartTime',Sign(1,1),'EndTime',Sign(end,1));
        end
        %------------------------------------------------------------------
        function Acc = resample1(Acc, nSamples)
            %Acc = RESAMPLE1(Acc, nSamples) resample the ACC object to a
            %'nSamples' samples time series, where 'nSamples' is a scalar.

            t=linspace(Acc.Time(1),Acc.Time(end),nSamples)';
            y=Acc.Data;
            while 1
                try
                    y=resample(y, nSamples, length(y));
                    break;
                catch
                    y=y(1:2:end);
                end
            end
            Name=Acc.Name;
            Acc=acc(t,y);
            Acc.Name=Name;
        end
        %------------------------------------------------------------------

        function SRS = srs(Acc, StaF, Q)
            %SRS = SRS(Acc, StaF, Q) plot the SRS curve of the time history,
            %where the 'StaF' is the starting frequency, and the 'Q' is the
            %quality factor.

            t=Acc.Time;
            y=Acc.Data;
            tmx=max(t);
            tmi=min(t);
            nnn = length(t);
            dt=(tmx-tmi)/(nnn-1);
            sr=1./dt;

            fn(1)=StaF;
            if fn(1)>sr/30.
                fn(1)=sr/30.;
            end
            damp=1./(2.*Q);
            j=1;
            while(1)
                if (fn(j) > sr/8.)
                    break
                end
                %fn(j+1)=fn(1)*(2. ^ (j*(1./1024.)));
                fn(j+1)=fn(1)*(2. ^ (j*(1./12.)));
                j=j+1;
            end
            %
            [a1,a2,b1,b2,b3]=srs.srs_coefficients(fn,damp,dt);

            tmax=(tmx-tmi) + 1./fn(1);

            limit = round( tmax/dt );
            tt=[t',tmx+dt:dt:tmax];
            yy=[y',zeros(size(y,2),limit-nnn)];

            ns=max(size(fn));
            %
            resp=zeros(length(fn),length(t));
            for j=1:ns
                %
                forward=[ b1(j),  b2(j),  b3(j) ];
                back   =[     1, -a1(j), -a2(j) ];
                %
                resp(j,:)=filter(forward,back,y',[],2);
                %
            end
            %maximaxSRS=max(x_pos,abs(x_neg));

            SRS=srs(fn, StaF, Q, t, resp);
            if nargout==0
                SRS.plot;
            end
        end

        function [maximaxSRS, fn] = srsm(Acc, StaF, Q)
            % SRS = srs(Acc, StaF, Q);
            t=Acc.Time;
            y=Acc.Data;
            tmx=max(t);
            tmi=min(t);
            nnn = length(t);
            dt=(tmx-tmi)/(nnn-1);
            sr=1./dt;

            fn(1)=StaF;
            if fn(1)>sr/30.
                fn(1)=sr/30.;
            end
            damp=1./(2.*Q);
            j=1;
            while(1)
                if (fn(j) > sr/8.)
                    break
                end
                %fn(j+1)=fn(1)*(2. ^ (j*(1./24.)));
                fn(j+1)=fn(1)*(2. ^ (j*(1./12.)));
                j=j+1;
            end
            %
            [a1,a2,b1,b2,b3]=srs.srs_coefficients(fn,damp,dt);

            ns=max(size(fn));
            %
            maximaxSRS=zeros(length(fn),size(y,2));

            maximaxSRS = gpuArray( single( maximaxSRS ));
            y = gpuArray( single(y) );
            for j=1:ns
                %
                forward=[ b1(j),  b2(j),  b3(j) ];
                back   =[     1, -a1(j), -a2(j) ];
                %
                z=filter(forward,back,y);
                maximaxSRS(j,:)=max(abs(z));
                %
            end
            maximaxSRS = gather( double( maximaxSRS ));
            fn = fn';
        end

        function plot(Acc)
            figure;
            plot(Acc.Time, Acc.Data, 'LineWidth', 0.5);

            xlabel({'Time (s)'});
            ylabel({'Acceleration (m/s^2)'});
            grid on;
            grid minor;
            set(gca,'FontSize', 14);
        end
        %------------------------------------------------------------------
        function [f, P1, P2, Y]=fft(Acc)
            L=length(Acc.Data);
            if rem(L,2)==0
                Y=fft(Acc.Data);
            else
                Y=fft(Acc.Data(1:end-1));
                L=L-1;
            end


            P2 = abs(Y/L);
            P1 = P2(1:L/2+1,:);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Acc.Sf*(0:(L/2))/L;
            if nargout==0
                figure;
                p1 = plot(f/1000,P1);
                p1.LineWidth = 1;
                xlabel('Frequency (kHz)');
                ylabel('Amplitude');
                set(gca,'FontSize', 12);
                grid on;
            end
        end
        %------------------------------------------------------------------
        function Acc=cumtrapz(Acc)
            Acc.Data=cumtrapz(Acc.Time,Acc.Data);
            if nargout==0
                Acc.plot;
            end
        end

        function Acc=diff(Acc)
            t=Acc.Time(1:end-1);
            y=diff(Acc.Data)*Acc.Sf;
            Acc=acc([t,y]);
            Acc.plot;
        end

        function Acc=extend(Acc0, duration)
            %Acc=EXTEND(Acc0, duration) extends the time history for
            %'duration' seconds.

            y=Acc0.Data;
            t=(0:1/Acc0.Sf:duration)';
            y=[y;zeros(length(t)-length(y), size(y,2))];
            Acc=acc(t,y);
        end

        function subplot(Acc, m, n)
            %SUBPLOT(Acc, m, n) first creates a m-by-n subplot, and then
            %plot all data columns in ACC subsequently.

            figure;
            for i=1:size(Acc.Data,2)
                subplot( m, n, i);
                plot(Acc.Time, Acc.Data(:,i));
            end
        end

    end

    methods (Static)

    end
end