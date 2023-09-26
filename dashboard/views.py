from django.contrib.auth.models import User
from django.views import generic
from blog.models import Posts, Comments
from blog.forms import SocialMediaModelForm
from dash import html, dcc
from django_plotly_dash import DjangoDash
import plotly.express as px
import pandas as pd

from django.contrib import messages
from django.shortcuts import redirect, get_object_or_404, reverse, render
from django.conf import settings
from accounts.models import Profile, SocialMedia, ContactModel
from accounts.forms import UserEditForm, ProfileEditForm, DeleteAccountForm
from django.http import HttpResponseRedirect
from django.contrib.auth.forms import PasswordChangeForm
from accounts.forms import ContactMessagesReplyForm, UserMessagesForm
from django.contrib.auth import update_session_auth_hash
from .models import Notifications
import dash_ag_grid as dag


class DashBoardUsers(generic.ListView):
    template_name = "dashboard/pages/all_users.html"
    model = User

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
        external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

        all_users_app = DjangoDash('all_users', external_stylesheets=external_stylesheets,
                                   external_scripts=external_scripts,
                                   add_bootstrap_links=True)

        all_users_df = pd.DataFrame({
            "Kullanıcı adı": [profile.username for profile in User.objects.all()],
            "Ad": [profile.first_name for profile in User.objects.all()],
            "Soyad": [profile.last_name for profile in User.objects.all()],
            "Email": [profile.email for profile in User.objects.all()],
            "Katılma Tarihi": [profile.date_joined.date() for profile in User.objects.all()],
            "Son oturum açma": [profile.last_login for profile in User.objects.all()],
        })

        all_users_app.layout = html.Div([

            html.Br(),

            html.H5("TÜM KULLANICILAR", className="text-success"),

            dag.AgGrid(
                id="all_users_app",
                style={'width': '100%'},
                rowData=all_users_df.to_dict("records"),
                columnSize="sizeToFit",
                defaultColDef={"resizable": True, "sortable": True, "filter": True, 'editable': True,
                               "minWidth": 125},
                dashGridOptions={
                    'pagination': True,
                    "rowSelection": "multiple",
                    "noRowsOverlayComponent": "CustomNoRowsOverlay",
                    "noRowsOverlayComponentParams": {
                        "message": "Gönderi bulunamadı",
                        "fontSize": 12,
                    },
                },
                columnDefs=[
                    {'field': 'Kullanıcı adı', 'headerName': 'Kullanıcı adı', 'filter': True},
                    {'field': 'Ad', 'headerName': 'Ad', 'filter': True},
                    {'field': 'Soyad', 'headerName': 'Soyad', 'filter': True},
                    {'field': 'Email', 'headerName': 'Email', 'filter': True},
                    {'field': 'Katılma Tarihi', 'headerName': 'Katılma Tarihi', 'filter': True},
                    {'field': 'Son oturum açma', 'headerName': 'Son oturum açma', 'filter': True},
                ]
            ),

            html.Div([
                dcc.Graph(figure=px.line(all_users_df.to_dict(), x="Katılma Tarihi", title="Kullanıcı",
                                         hover_data=["Kullanıcı adı", "Ad", 'Email'], markers=True))

            ])

        ], style={"margin": 20}, )

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


class DashBoardPosts(generic.ListView):
    template_name = "dashboard/pages/all_posts.html"
    model = Posts

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        external_stylesheets = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css']
        external_scripts = ['https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.bundle.min.js']

        all_post_app = DjangoDash('all_posts', external_stylesheets=external_stylesheets,
                                  external_scripts=external_scripts,
                                  add_bootstrap_links=True)

        post_df = pd.DataFrame({
            "#": [post.id for post in Posts.objects.all()],
            "Tarih": [post.created.date() for post in Posts.objects.all()],
            "Başlık": [post.title[:20] for post in Posts.objects.all()],
            "Yazar": [post.author.username for post in Posts.objects.all()],
            "Kategori": [post.category.title for post in Posts.objects.all()],
            "Likes": [post.likes.count() for post in Posts.objects.all()],
            "Dislikes": [post.dislike.count() for post in Posts.objects.all()],
        })

        all_post_app.layout = html.Div([

            html.Br(),

            html.H5("TÜM GÖNDERİLER", className="text-success"),

            dag.AgGrid(
                id="frame_seq_table",
                style={'width': '100%'},
                rowData=post_df.to_dict("records"),
                columnSize="sizeToFit",
                defaultColDef={"resizable": True, "sortable": True, "filter": True, 'editable': True,
                               "minWidth": 125},
                dashGridOptions={
                    'pagination': True,
                    "rowSelection": "multiple",
                    "noRowsOverlayComponent": "CustomNoRowsOverlay",
                    "noRowsOverlayComponentParams": {
                        "message": "Gönderi bulunamadı",
                        "fontSize": 12,
                    },
                },
                columnDefs=[
                    {'field': '#', 'headerName': 'İD', 'filter': True},
                    {'field': 'Tarih', 'headerName': 'Tarih', 'filter': True},
                    {'field': 'Başlık', 'headerName': 'Başlık', 'filter': True, },
                    {'field': 'Yazar', 'headerName': 'Yazar', 'filter': True, },
                    {'field': 'Kategori', 'headerName': 'Kategori', 'filter': True},
                    {'field': 'Likes', 'headerName': 'Likes', 'filter': True},
                    {'field': 'Dislikes', 'headerName': 'Dislikes', 'filter': True},
                ]
            ),

            html.Div([
                dcc.Graph(figure=px.bar(data_frame=post_df.to_dict(), x="Kategori", color="Kategori",
                                        title="Gönderi-Kategori", labels={'count': 'Gönderi Sayısı'}, barmode="group",
                                        hover_data=['Tarih', 'Başlık', 'Likes'], text_auto=True, ))
            ])

        ], style={"margin": 20}, )

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


def profile_update(request, pk, username):
    if not request.user.username == username:
        messages.error(request, "Yetkili girişi yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

    user = Profile.objects.get(user=pk)
    user_edit_form = UserEditForm(request.POST or None, instance=request.user)
    profile_edit_form = ProfileEditForm(request.POST or None, request.FILES or None, instance=request.user.profile)
    password_form = PasswordChangeForm(request.user, data=request.POST or None)
    delete_account_form = DeleteAccountForm(request.POST or None)
    notifications = Notifications.objects.filter(user=request.user, is_read=False).order_by("-created")
    messages_list = ContactModel.objects.filter(receiver=request.user.pk)
    messages_count = ContactModel.objects.filter(receiver=request.user.pk, is_read=False).count()

    if request.method == "POST":

        if user_edit_form.is_valid():

            user_edit_form.save()
            messages.success(request, 'Kullanıcı bilgileriniz başarılı bir şekilde güncellendi')
            return redirect(request.META['HTTP_REFERER'])

        elif profile_edit_form.is_valid():
            profile_edit_form.save()
            messages.success(request, 'Profil bilgileriniz başarılı bir şekilde güncellendi')
            return redirect(request.META['HTTP_REFERER'])

        elif password_form.is_valid():

            update_session_auth_hash(request, password_form.user)
            password_form.save()
            messages.success(request, 'Kullanıcı şifreniz başarılı bir şekilde güncellendi')
            return redirect(request.META['HTTP_REFERER'])

        elif delete_account_form.is_valid():
            confirm = delete_account_form.cleaned_data['confirm']
            if confirm is True:
                user.delete()
                messages.success(request, 'Profil bilgileriniz başarılı silinmiştir')
                return redirect('login')
        else:
            profile_edit_form = ProfileEditForm(instance=request.user.profile)
            user_edit_form = UserEditForm(instance=request.user)

    return render(request, 'dashboard/pages/update.html',
                  {'user_edit_form': user_edit_form, 'profile_edit_form': profile_edit_form,
                   'password_form': password_form, 'delete_account_form': delete_account_form,
                   'messages_list': messages_list, 'messages_count': messages_count, 'notifications': notifications})


class MessagesListView(generic.ListView):
    template_name = 'dashboard/pages/messages.html'
    model = ContactModel
    context_object_name = 'messages_list'

    def get_queryset(self):
        return ContactModel.objects.filter(receiver=self.request.user)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


class MessagesDetailView(generic.DetailView):
    template_name = 'dashboard/pages/messages.html'
    model = ContactModel

    def get(self, request, *args, **kwargs):
        ContactModel.objects.update(is_read=True)
        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


def mark_read_message(request, pk):
    message = get_object_or_404(ContactModel, pk=pk)
    message.is_read = True
    message.save()
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def mark_as_read_all(request):
    for message in ContactModel.objects.filter(receiver=request.user.pk).all():
        message.is_read = True
        message.save()
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def report_message(request, pk):
    message = get_object_or_404(ContactModel, pk=pk)
    message.is_report = True
    message.save()
    messages.success(request, "Mesajı rapor ettiniz!")
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


def message_delete(request, pk):
    obj = get_object_or_404(ContactModel, pk=pk)
    obj.delete()
    messages.success(request, "Mesaj başarılı bir şekilde silindi!")
    return HttpResponseRedirect(reverse('dashboard:messages_list', kwargs={
        'user': request.user.id
    }))


class SocialMediaListView(generic.ListView):
    template_name = 'dashboard/pages/index.html'
    model = SocialMedia

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.pk)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


class SocialMediaCreateView(generic.CreateView):
    template_name = 'dashboard/pages/create.html'
    model = SocialMedia
    form_class = SocialMediaModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.id)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user_id = self.request.user.id
        if SocialMedia.objects.filter(user=self.request.user.id, social=instance.social).exists():
            messages.error(self.request, "Sosyal Medya bilgisi mevcut!")
            return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
                'user': self.request.user
            }))
        instance.save()
        messages.success(self.request, "Sosyal Medya Bilgisi başarılı bir şekilde Eklendi!")
        return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
            'user': self.request.user
        }))

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


class SocialMediaUpdateView(generic.UpdateView):
    template_name = 'dashboard/pages/update.html'
    model = SocialMedia
    form_class = SocialMediaModelForm

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return SocialMedia.objects.filter(user=self.request.user.id)

    def form_valid(self, form):
        instance = form.save(commit=False)
        instance.user_id = self.request.user.id
        instance.save()
        messages.success(self.request, "Sosyal Medya Bilgisi başarılı bir şekilde güncellendi!")
        return HttpResponseRedirect(reverse('dashboard:social_dash', kwargs={
            'user': self.request.user
        }))

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


def user_reply_message(request, pk, username, user_pk):
    if not request.user.is_staff:
        messages.error(request, "Yetkili girişi yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

    form = ContactMessagesReplyForm(request.POST or None)
    contact_object = ContactModel.objects.get(pk=pk)
    receiver = User.objects.get(pk=user_pk)

    notifications = Notifications.objects.filter(user=request.user, is_read=False).order_by("-created")
    messages_list = ContactModel.objects.filter(receiver=request.user.pk)
    messages_count = ContactModel.objects.filter(receiver=request.user.pk, is_read=False).count()

    if request.method == "POST":
        if form.is_valid():
            sender = request.user
            message = form.cleaned_data['content']
            ContactModel.objects.create(sender=sender, receiver=receiver, title=contact_object.title, content=message,
                                        contact_email=receiver.email)
            messages.success(request, 'Mesajnız gönderilmiştir...')

            return HttpResponseRedirect(reverse('dashboard:messages_list',
                                                kwargs={'user': request.user.username}))
        else:
            form = ContactMessagesReplyForm()

    return render(request, 'dashboard/pages/messages.html', {
        'form': form, 'contact_object': contact_object,
        'messages_list': messages_list,
        'messages_count': messages_count,
        'notifications': notifications
    })


def user_sent_message(request, pk, username):
    if not request.user.is_staff:
        messages.error(request, "Yetkili girişi yapınız!")
        return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
    form = UserMessagesForm(request.POST or None)

    notifications = Notifications.objects.filter(user=request.user, is_read=False).order_by("-created")
    messages_list = ContactModel.objects.filter(receiver=request.user.pk)
    messages_count = ContactModel.objects.filter(receiver=request.user.pk, is_read=False).count()

    if request.method == "POST":
        if form.is_valid():

            title = form.cleaned_data['title']
            message = form.cleaned_data['content']
            contact_email = form.cleaned_data['contact_email']
            sender = request.user
            receiver = User.objects.get(username=form.cleaned_data['receiver'])

            if str(request.user.username) == str(form.cleaned_data['receiver']):
                messages.error(request, 'Kendinize mesaj gönderdiniz. Farklı bir kullanıcı seçiniz.')
                return redirect(request.META['HTTP_REFERER'])

            ContactModel.objects.create(title=title, sender=sender, receiver=receiver, content=message,
                                        contact_email=contact_email)
            messages.success(request, 'Mesajnız gönderilmiştir...')
            return HttpResponseRedirect(reverse('dashboard:messages_list',
                                                kwargs={'user': request.user}))

        else:
            form = UserMessagesForm()

    return render(request, 'dashboard/pages/messages.html',
                  {
                      'form': form,
                      'messages_list': messages_list,
                      'messages_count': messages_count,
                      'notifications': notifications
                  })


class MyPostList(generic.ListView):
    template_name = "dashboard/pages/index.html"
    model = Posts
    context_object_name = "my_posts"

    def get(self, request, *args, **kwargs):
        if not self.request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))
        return super().get(request, *args, **kwargs)

    def get_queryset(self):
        return Posts.objects.filter(author=self.request.user).order_by('-created')

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context


class NotificationsListView(generic.ListView):
    template_name = "dashboard/pages/notifications.html"
    model = Notifications

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['notifications'] = Notifications.objects.filter(user=self.request.user, is_read=False).order_by(
            "-created")
        context['messages_list'] = ContactModel.objects.filter(receiver=self.request.user.pk)
        context['messages_count'] = ContactModel.objects.filter(receiver=self.request.user.pk, is_read=False).count()
        return context
